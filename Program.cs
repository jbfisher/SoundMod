using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Numerics;
using System.Media;
using System.Runtime.InteropServices;
using System.Threading;

namespace DemodConsole
{
	class Program
	{
		private bool Running = true;

		const uint UPPER_PRIME_LIMIT = (44100 / 2);
		const uint FREQ_BOUND_LOW = 3500;
		const uint FREQ_BOUND_HIGH = 8000;
		const uint NUM_FREQS = 39;
		const uint FREQS_PER = 11; // (NUM_FREQS - 6) / NUM_VALS
		const uint NUM_VALS = 3;

		uint START_NUM = 20000;
		uint RecordingNumber;

		int AUDIO_FILE_DURATION = 200;
		int DEMODULATOR_PAUSE = 160;

		double StdDevFactor = 0.8;

		string FileName = "Rec";
		string FileExt = ".wav";

		ushort FrequencySpacing;
		ushort AdjustedFrequencySpacing;
		ushort[] TargetFrequency;
		ushort[] AdjustedTargetFrequency;

		[DllImport("winmm.dll", EntryPoint = "mciSendStringA", CharSet = CharSet.Ansi, SetLastError = true, ExactSpelling = true)]
		private static extern int mciSendString(string lpstrCommand, string lpstrReturnString, int uReturnLength, int hwndCallback);

		static void Main(string[] args)
		{
			Program P = new Program();
			P.InitializeFrequencies();
			P.DeletePreviousAudioFiles();

			ThreadStart TS1 = new ThreadStart(P.ListenAndRecord);
			Thread ListenerThread = new Thread(TS1);
			ListenerThread.Start();

			ThreadStart TS2 = new ThreadStart(P.Demodulate);
			Thread DemodulatorThread = new Thread(TS2);
			DemodulatorThread.Start();

			ListenerThread.Join();
			DemodulatorThread.Join();
		}

		public void ListenAndRecord()
		{
			RecordingNumber = START_NUM;
			int Error;
			while (Running)
			{
				string FullFileName = FileName + RecordingNumber.ToString() + FileExt;
				Error = mciSendString("open new type waveaudio alias recsound", "", 0, 0);
				if (Error != 0) Console.Out.WriteLine("mciSendString open failed; error: " + Error.ToString());

				Error = mciSendString("set recsound bitspersample 16 channels 1 samplespersec 44100 alignment 2 bytespersec 88200", "", 0, 0);
				if (Error != 0) Console.Out.WriteLine("mciSendString set failed; error: " + Error.ToString());

				Error = mciSendString("record recsound", "", 0, 0);
				if (Error != 0) Console.Out.WriteLine("mciSendString recsound failed; error: " + Error.ToString());
				Thread.Sleep(AUDIO_FILE_DURATION);
				Error = mciSendString(@"save recsound " + FullFileName, "", 0, 0);
				if (Error != 0) Console.Out.WriteLine("mciSendString save failed; error: " + Error.ToString());
				Error = mciSendString("close recsound ", "", 0, 0);
				if (Error != 0) Console.Out.WriteLine("mciSendString close failed; error: " + Error.ToString());
				//else
					//Console.Out.WriteLine("Saved a WAV file.");
				RecordingNumber++;
			}
		}

		private void InitializeFrequencies()
		{
			// Calculate all primes up to a limit; true values in the array will indicate primeness
			bool[] Primes = new bool[UPPER_PRIME_LIMIT];

			for (int PossPrime = 0; PossPrime < UPPER_PRIME_LIMIT; PossPrime++)
				Primes[PossPrime] = true;

			Primes[0] = false;
			Primes[1] = false;

			// Remove all multiples; remaining true values will be prime
			for (int PossPrime = 0; PossPrime < UPPER_PRIME_LIMIT / 2; PossPrime++)
			{
				if (Primes[PossPrime])
				{
					for (int Multiple = PossPrime * 2; Multiple < UPPER_PRIME_LIMIT; Multiple += PossPrime)
						Primes[Multiple] = false;
				}
			}

			TargetFrequency = new ushort[NUM_FREQS];
			AdjustedTargetFrequency = new ushort[NUM_FREQS];

			// We will take frequency values that almost spaced evenly within the specified range
			// between FREQ_BOUND_HIGH and FREQ_BOUND_LOW, but jump up to the next prime value
			FrequencySpacing = (ushort)((FREQ_BOUND_HIGH - FREQ_BOUND_LOW) / NUM_FREQS);

			ushort NextFrequency = (ushort)FREQ_BOUND_LOW;
			ushort FrequencyToCheck = NextFrequency;
			for (uint i = 0; i < NUM_FREQS; i++)
			{
				while (!Primes[FrequencyToCheck])
					FrequencyToCheck++;
				TargetFrequency[i] = FrequencyToCheck;

				NextFrequency += FrequencySpacing;
				FrequencyToCheck = NextFrequency;

				Console.Out.WriteLine(TargetFrequency[i].ToString());
			}
		}

		private void DeletePreviousAudioFiles()
		{
			uint RecordingNumber = START_NUM;
			string ThisFileName = FileName + RecordingNumber.ToString() + FileExt;
			while (File.Exists(ThisFileName))
			{
				File.Delete(ThisFileName);
				RecordingNumber++;
				ThisFileName = FileName + RecordingNumber.ToString() + FileExt;
			}
		}

		private void Demodulate()
		{
			Console.Out.WriteLine("Beginning demodulation.");
			uint FileNumberToRead;
			uint LastFileRead = 0;
			string ThisFileName;
			FileStream FS;
			BinaryReader BR;
			short LastSignalReceived = -1;
			char[][][] DemoddedValues = new char[3][][];
			for (uint i = 0; i < 3; i++)
			{
				DemoddedValues[i] = new char[NUM_VALS][];
				for (uint ii = 0; ii < NUM_VALS; ii++)
					DemoddedValues[i][ii] = new char[FREQS_PER];
			}
			while (Running)
			{
				FileNumberToRead = RecordingNumber;
				FileNumberToRead--;
				ThisFileName = FileName + FileNumberToRead.ToString() + FileExt;
				// Read file
				//ThisFileName = "C:\\Users\\J. Fisher\\Desktop\\AModTest.wav";
				//ThisFileName = "C:\\Users\\J. Fisher\\Desktop\\Rec10006.wav";
				//Running = false;

				if ((FileNumberToRead >= START_NUM && FileNumberToRead != LastFileRead) || !Running)
				//if (NumberToRead >= START_NUM || !Running)
					{
					FS = new FileStream(ThisFileName, FileMode.Open, FileAccess.Read);
					BR = new BinaryReader(FS);
					//Console.Out.WriteLine();
					//Console.Out.WriteLine("Reading file: " + ThisFileName);
					if (!Running)
						Console.Out.WriteLine();  

					uint FileSize;

					uint FormatChunkSize;
					ushort FormatTag;
					ushort NumChannels;
					uint SamplesPerSecond;
					uint BytesPerSecond;
					ushort BlockAlign;
					ushort BitsPerSample;

					uint DataSize;
					uint NumberOfSamples;
					float Duration;
					float[] Data1;
					float[] Data2 = null;
					short TempData;

					// Read header
					for (uint i = 0; i < 4; i++)
						BR.ReadChar().ToString();
					FileSize = BR.ReadUInt32();
					for (uint i = 0; i < 4; i++)
						BR.ReadChar().ToString();

					// Read format chunk
					for (uint i = 0; i < 4; i++)
						BR.ReadChar().ToString();
					FormatChunkSize = BR.ReadUInt32();
					FormatTag = BR.ReadUInt16();
					NumChannels = BR.ReadUInt16();
					SamplesPerSecond = BR.ReadUInt32();
					BytesPerSecond = BR.ReadUInt32();
					BlockAlign = BR.ReadUInt16();
					BitsPerSample = BR.ReadUInt16();

					// Skip any other chunks before data, then read data
					char TempChar = ' ';
					while (TempChar != 'd')
						TempChar = BR.ReadChar();
					for (uint i = 0; i < 3; i++)
						BR.ReadChar().ToString();
					DataSize = BR.ReadUInt32();
					NumberOfSamples = DataSize / BlockAlign;
					Data1 = new float[NumberOfSamples];
					if (NumChannels == 2)
						Data2 = new float[NumberOfSamples];
					Duration = (float)NumberOfSamples / (float)SamplesPerSecond;

/*					// Read header
					for (uint i = 0; i < 4; i++)
						Console.Out.Write(BR.ReadChar().ToString());
					Console.Out.Write("\t");
					FileSize = BR.ReadUInt32();
					Console.Out.Write(FileSize.ToString() + "\t");
					for (uint i = 0; i < 4; i++)
						Console.Out.Write(BR.ReadChar().ToString());
					Console.Out.Write("\t");

					// Read format chunk
					for (uint i = 0; i < 4; i++)
						Console.Out.Write(BR.ReadChar().ToString());
					Console.Out.WriteLine();
					FormatChunkSize = BR.ReadUInt32();

					Console.Out.Write("FormatChunkSize: " + FormatChunkSize.ToString() + "\t");
					FormatTag = BR.ReadUInt16();
					Console.Out.Write("FormatTag: " + FormatTag.ToString() + "\t");
					NumChannels = BR.ReadUInt16();
					Console.Out.Write("NumChannels: " + NumChannels.ToString() + "\t");
					SamplesPerSecond = BR.ReadUInt32();
					Console.Out.Write("SamplesPerSecond: " + SamplesPerSecond.ToString() + "\t");
					BytesPerSecond = BR.ReadUInt32();
					Console.Out.Write("BytesPerSecond: " + BytesPerSecond.ToString() + "\t");
					BlockAlign = BR.ReadUInt16();
					Console.Out.Write("BlockAlign: " + BlockAlign.ToString() + "\t");
					BitsPerSample = BR.ReadUInt16();
					Console.Out.Write("BitsPerSample: " + BitsPerSample.ToString());
					Console.Out.WriteLine();

					// Skip any other chunks before data, then read data
					char TempChar = ' ';
					while (TempChar != 'd')
						TempChar = BR.ReadChar();
					Console.Out.Write(TempChar.ToString());
					for (uint i = 0; i < 3; i++)
						Console.Out.Write(BR.ReadChar().ToString());
					Console.Out.Write("\t");
					DataSize = BR.ReadUInt32();
					Console.Out.Write("Size: " + DataSize.ToString() + "\t");
					NumberOfSamples = DataSize / BlockAlign;
					Data1 = new float[NumberOfSamples];
					if (NumChannels == 2)
						Data2 = new float[NumberOfSamples];
					Console.Out.Write("# of samples*: " + NumberOfSamples.ToString() + "\t");
					Duration = (float)NumberOfSamples / (float)SamplesPerSecond;
					Console.Out.Write("Duration*: " + Duration.ToString() + " sec");
					Console.Out.Write("\t");
*/
					for (uint i = 0; i < NumberOfSamples; i++)
					{
						TempData = BR.ReadInt16();
						Data1[i] = (float)TempData;
						if (NumChannels == 2)
						{
							TempData = BR.ReadInt16();
							Data2[i] = (float)TempData;
							Console.Out.Write("`");
						}
					}
					//Console.Out.WriteLine("Done!");
					FS.Close();

					// Pack extra info to make sure the number of samples is enough for the algorithm.
					// ??? May be necessary if the sample array is too short
					// Not implemented right now

					Complex[] FFTVal = new Complex[NumberOfSamples];				// storage for FFT answer

					for (uint i = 0; i < NumberOfSamples; i++)
						FFTVal[i] = new Complex(Data1[i], 0.0);		// copy into FFTVal[] for FFT work & result

					FFT2(FFTVal, NumberOfSamples);

					uint NumberOfFFTVals = NumberOfSamples / 2;
					// Convert values to simple float type
					float[] FinalValues = new float[NumberOfFFTVals];
					for (uint i = 0; i < NumberOfFFTVals; i++)
						FinalValues[i] = (float)Complex.Abs(FFTVal[i]);

					float FreqResolution = SamplesPerSecond / (float)NumberOfSamples; // freq step in FFT result
/*
					// Write to file for later inspection
					FileStream FOut = new FileStream("FFTVal0.txt", FileMode.Create);
					StreamWriter SOut = new StreamWriter(FOut);
					float X;
					for (uint i = 0; i < NumberOfFFTVals; i++)
					{
						X = i * FreqResolution;
						SOut.WriteLine(X.ToString() + "\t" + FinalValues[i]);
					}
*/ 
					// Calculate standard deviation
					double FFTMean = 0.0;
					for (uint i = 0; i < NumberOfFFTVals; i++)
						FFTMean += FinalValues[i];
					FFTMean /= NumberOfFFTVals;
					double FFTSum = 0.0;
					for (uint i = 0; i < NumberOfFFTVals; i++)
						FFTSum += Math.Pow(FinalValues[i] - FFTMean, 2.0);
					FFTSum /= (NumberOfFFTVals - 1);
					double StdDev = Math.Sqrt(FFTSum);

					// Shift values according to FreqResolution
					AdjustedFrequencySpacing = (ushort)((float)FrequencySpacing / FreqResolution);
					for (uint i = 0; i < NUM_FREQS; i++)
						AdjustedTargetFrequency[i] = (ushort)(((float)TargetFrequency[i]) / FreqResolution);

					ushort CheckBandWidth = (ushort)(AdjustedFrequencySpacing / 2);

					ushort UpperLimit, LowerLimit;
					double Amplitude;

					short TransmissionPiece = -1;

					// Check for second-highest frequency
					ushort FirstTransmissionBit = 0;
					for (ushort Bit = 0; Bit < 3; Bit++)
					{
						LowerLimit = (ushort)(AdjustedTargetFrequency[NUM_FREQS - 6 + Bit] - (CheckBandWidth / 2));
						UpperLimit = (ushort)(AdjustedTargetFrequency[NUM_FREQS - 6 + Bit] + (CheckBandWidth / 2));
						Amplitude = 0.0;
						for (int i = LowerLimit; i < UpperLimit; i++)
							Amplitude += FinalValues[i];
						Amplitude /= CheckBandWidth;
						if (Amplitude > StdDev * StdDevFactor)
							FirstTransmissionBit++;
					}
					if (FirstTransmissionBit > 1)
						TransmissionPiece += 1;

					// Check for highest frequency
					ushort SecondTransmissionBit = 0;
					for (ushort Bit = 0; Bit < 3; Bit++)
					{
						LowerLimit = (ushort)(AdjustedTargetFrequency[NUM_FREQS - 3 + Bit] - (CheckBandWidth / 2));
						UpperLimit = (ushort)(AdjustedTargetFrequency[NUM_FREQS - 3 + Bit] + (CheckBandWidth / 2));
						Amplitude = 0.0;
						for (int i = LowerLimit; i < UpperLimit; i++)
							Amplitude += FinalValues[i];
						Amplitude /= CheckBandWidth;
						if (Amplitude > StdDev * StdDevFactor)
							SecondTransmissionBit++;
					}
					if (SecondTransmissionBit > 1)
						TransmissionPiece += 2;

					if (TransmissionPiece >= 0)
					{
//						Console.Out.WriteLine(); 
//						Console.Out.Write("Signal: " + TransmissionPiece.ToString());
						string TempString;
						short Value;
						for (uint Val = 0; Val < NUM_VALS; Val++)
						{
							Value = 0;
							TempString = "";
							for (uint Freq = 0; Freq < FREQS_PER; Freq++)
							{
								// Check to see if frequency is present
								LowerLimit = (ushort)(AdjustedTargetFrequency[Val * FREQS_PER + Freq] - (CheckBandWidth / 2));
								UpperLimit = (ushort)(AdjustedTargetFrequency[Val * FREQS_PER + Freq] + (CheckBandWidth / 2));
								Amplitude = 0.0;
								for (int i = LowerLimit; i < UpperLimit; i++)
									Amplitude += FinalValues[i];
								Amplitude /= CheckBandWidth;
								//if (Amplitude > BadThreshold)
								if (Amplitude > StdDev * StdDevFactor)
								// Frequency is present
								{
									TempString += "1";
									if (Freq != (FREQS_PER - 1))
										Value += (short)Math.Pow(2, Freq);
									else
										Value *= -1;
								}
								else
								{
									// Frequency is not present
									TempString += "0";
								}
							}
							for (int BitPosition = 0; BitPosition < FREQS_PER; BitPosition++)
								DemoddedValues[TransmissionPiece][(Val + TransmissionPiece) % NUM_VALS][BitPosition] = TempString[BitPosition];
//							Console.Out.Write(TempString + " " + Value.ToString() + "\t");
						}
						Console.Out.WriteLine();
						if (TransmissionPiece == 0)
							LastSignalReceived = 0;
						else if (TransmissionPiece == 1)
						{
							if (LastSignalReceived == 0)
								LastSignalReceived = 1;
							else
								LastSignalReceived = -1;
						}
						else if (TransmissionPiece == 2)
						{
							if (LastSignalReceived == 1)
							{
								// Display to screen
								for (uint Transmission = 0; Transmission < 3; Transmission++)
								{
									for (uint Val = 0; Val < NUM_VALS; Val++)
									{
										for (uint Freq = 0; Freq < FREQS_PER; Freq++)
											Console.Out.Write(DemoddedValues[Transmission][Val][Freq]);
										Console.Out.Write(" ");
									}
									Console.Out.WriteLine();
								}

								// Calculate bit-corrected values
								for (uint Val = 0; Val < NUM_VALS; Val++)
								{
									short ValueTotal = 0;
									for (uint Freq = 0; Freq < FREQS_PER; Freq++)
									{
										ushort BitSum = 0;
										for (uint Transmission = 0; Transmission < 3; Transmission++)
											if (DemoddedValues[Transmission][Val][Freq] == '1')
												BitSum++;
										// Rewrite over first set of values
										if (BitSum > 1)
										{
											if (Freq != FREQS_PER - 1)
												ValueTotal += (short)Math.Pow(2, Freq);
											else
												ValueTotal *= -1;
										}
									}
									Console.Out.Write(ValueTotal.ToString() + "\t");
								}
								Console.Out.WriteLine();
							}
							LastSignalReceived = -1;
						}
					}
					// Delete file
					//					File.Delete(ThisFileName);
					LastFileRead = FileNumberToRead;
					if (!Running)
						Console.Out.WriteLine();
				}
				else
				{
					//MessageBox.Show("File not present: " + ThisFileName);
				}
				// Wait
				Console.Out.Write(".");
				Thread.Sleep(DEMODULATOR_PAUSE);
			}
		}

		// Following code adapted from C code on wikipedia.org under Fast Fourier Transform:
		// https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

		/* fft.cpp
		 * 
		 * This is a KISS implementation of
		 * the Cooley-Tukey recursive FFT algorithm.
		 * This works, and is visibly clear about what is happening where.
		 */

		// Separate even/odd elements to lower/upper halves of array respectively.
		// Due to Butterfly combinations, this turns out to be the simplest way 
		// to get the job done without clobbering the wrong elements.
		private void Separate(Complex[] A, uint N)
		{
			Complex[] B = new Complex[N / 2];  // get temp heap storage
			for (uint i = 0; i < N / 2; i++)	// copy all odd elements to heap storage
				B[i] = A[i * 2 + 1];
			for (uint i = 0; i < N / 2; i++)	// copy all even elements to lower-half of a[]
				A[i] = A[i * 2];
			for (uint i = 0; i < N / 2; i++)	// copy all odd (from heap) to upper-half of a[]
				A[i + N / 2] = B[i];
		}

		// N must be a power-of-2, or bad things will happen.
		// Currently no check for this condition.
		//
		// N input samples in X[] are FFT'd and results left in X[].
		// Because of Nyquist theorem, N samples means 
		// only first N/2 FFT results in X[] are the answer.
		// (upper half of X[] is a reflection with no new information).
		private void FFT2(Complex[] X, uint N)
		{
			if (N < 2) { }	// Bottom of recursion; do nothing here, because already X[0] = x[0]
			else
			{
				Separate(X, N);	  // all evens to lower half, all odds to upper half
				FFT2(X, N / 2);   // recurse even items
				Complex[] TempX = new Complex[N / 2];
				for (uint i = 0; i < N / 2; i++)
					TempX[i] = X[i + N / 2];
				FFT2(TempX, N / 2);   // recurse odd  items
				for (uint i = 0; i < N / 2; i++)
					X[i + N / 2] = TempX[i];
				// combine results of two half recursions
				for (uint k = 0; k < N / 2; k++)
				{
					Complex e = X[k];   // even
					Complex o = X[k + N / 2];   // odd
					// w is the "twiddle-factor"
					Complex w = new Complex(Math.Cos(-2.0 * Math.PI * k / N), Math.Sin(-2.0 * Math.PI * k / N));
					X[k] = e + w * o;
					X[k + N / 2] = e - w * o;
				}
			}
		}
	}
}

