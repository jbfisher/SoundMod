// AudMod.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <time.h>

using namespace std;

const int M_PI = 3.14159265359;



//#define NUM_FREQS 11

const int UPPER_PRIME_LIMIT = (44100 / 2);
const int FREQ_BOUND_LOW = 18000;
const int FREQ_BOUND_HIGH = 21000;
#define NUM_FREQS 45
#define FREQS_PER 11 // (NUM_FREQS - 1) / 4

int _tmain(int argc, _TCHAR* argv[])
{
	srand(time(NULL));
	const int BITS_PER_BYTE = 8;



	// Calculate all primes up to a limit; true values in the array will indicate primeness
	bool Primes[UPPER_PRIME_LIMIT];

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

	unsigned short Frequency[NUM_FREQS];

	// We will take frequency values that almost spaced evenly within the specified range
	// between FREQ_BOUND_HIGH and FREQ_BOUND_LOW, but jump up to the next prime value
	int Spacing = (FREQ_BOUND_HIGH - FREQ_BOUND_LOW) / NUM_FREQS;

	int NextFrequency = FREQ_BOUND_LOW;
	int FrequencyToCheck = NextFrequency;
	for (int i = 0; i < NUM_FREQS; i++)
	{
		while (!Primes[FrequencyToCheck])
			FrequencyToCheck++;
		Frequency[i] = FrequencyToCheck;

		NextFrequency += Spacing;
		FrequencyToCheck = NextFrequency;

		cout << Frequency[i] << endl;
	}



	ofstream OutFile("AModTest.wav", ios::binary | ios::out);
	
	unsigned int FileSize;

	unsigned int FormatChunkSize = 16;
	unsigned short FormatTag = 1;
	unsigned short NumChannels = 1;
	unsigned int SamplesPerSecond = 44100;
	unsigned short BytesPerSample = sizeof(short);
	unsigned int BytesPerSecond = SamplesPerSecond * BytesPerSample;
	unsigned short BlockAlign = NumChannels * BytesPerSample;
	unsigned short BitsPerSample = BytesPerSample * BITS_PER_BYTE;

	float Duration = 0.5;
	unsigned int DataSize = (int)(BytesPerSecond * Duration);

	FileSize = 8 + 8 + FormatChunkSize + 8 + DataSize;



	float Wavelength[NUM_FREQS];
	for (int i = 0; i < NUM_FREQS; i++)
		Wavelength[i] = SamplesPerSecond / (float)Frequency[i];

	float Noise = 0.25;

	short NextData;
	short MaxAmplitude = 32760;
	float BaseValue;


	short TimeStamp = 5;
	short X = -105, Y = 40, Z = 1000;
//	short TimeStamp = -1023, X = -1023, Y = -1023, Z = -1023; // all 1s in the signal

	short Values[4];
	Values[0] = TimeStamp;
	Values[1] = X;
	Values[2] = Y;
	Values[3] = Z;
	for (int NextVal = 0; NextVal < 4; NextVal++)
		cout << "Value: " << Values[NextVal] << endl;

	char Bits[4][FREQS_PER];

	for (int NextVal = 0; NextVal < 4; NextVal++)
	{
		// Default to 0s
		for (int i = 0; i < FREQS_PER; i++)
			Bits[NextVal][i] = '0';

		// Encode a negative value with a trailing 1
		if (Values[NextVal] < 0)
		{
			Bits[NextVal][FREQS_PER - 1] = '1';
			Values[NextVal] = Values[NextVal] * -1;
		}
		// Encode values in binary as characters
		for (int Place = 0; Place < (FREQS_PER - 1); Place++)
		{
			if (Values[NextVal] % 2 == 1)
				Bits[NextVal][Place] = '1';
			Values[NextVal] = Values[NextVal] / 2;
		}
	}
	for (int NextVal = 0; NextVal < 4; NextVal++)
		cout << Bits[NextVal] << endl;

	// Write header
	OutFile.write("RIFF", 4);
	OutFile.write((char*)&FileSize, sizeof(FileSize));
	OutFile.write("WAVE", 4);

	// Write format chunk
	OutFile.write("fmt ", 4);
	OutFile.write((char*)&FormatChunkSize, sizeof(FormatChunkSize));
	OutFile.write((char*)&FormatTag, sizeof(FormatTag));
	OutFile.write((char*)&NumChannels, sizeof(NumChannels));
	OutFile.write((char*)&SamplesPerSecond, sizeof(SamplesPerSecond));
	OutFile.write((char*)&BytesPerSecond, sizeof(BytesPerSecond));
	OutFile.write((char*)&BlockAlign, sizeof(BlockAlign));
	OutFile.write((char*)&BitsPerSample, sizeof(BitsPerSample));

	// Write data chunk
	OutFile.write("data", 4);
	OutFile.write((char*)&DataSize, sizeof(DataSize));
	int UpperLimit = DataSize / 2;
	for (int Sample = 0; Sample < UpperLimit; Sample++)
	{
		BaseValue = sin(Sample / Wavelength[NUM_FREQS - 1] * 2 * M_PI);
		for (int NextVal = 0; NextVal < 4; NextVal++)
			for (int i = 0; i < FREQS_PER; i++)
				if (Bits[NextVal][i] == '1')
					BaseValue += sin(Sample / Wavelength[NextVal * FREQS_PER + i] * 2 * M_PI);
		BaseValue = BaseValue / (double)NUM_FREQS;
		NextData = (short)(BaseValue * MaxAmplitude);
		NextData += (rand() % (32670 * 2) - 32670) * Noise;
		OutFile.write((char*)&NextData, sizeof(NextData));
	}
	OutFile.close();

	return 0;
}

