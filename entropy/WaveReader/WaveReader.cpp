// WaveReader.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include "windows.h"
#include "mmsystem.h"
#include <fstream>
#include <time.h>

using namespace std;

struct wave_header
{
	char chunkID[4]; //ID RIFF
	unsigned long chunkSize; //rozmiar danych
	char format[4]; //format pliku (WAV)
	char subchunk1ID[4];  //pocz¹tek czêœci opisowej (fmt)
	unsigned long subchunk1Size; //rozmiar czêœci opisowej, dla fmt zwykle 16
	unsigned short audioFormat; //rodzaj kompresji, 1 - bez kompresji, modulacja PCM
	unsigned short numChannels; //liczba kana³ów
	unsigned long sampleRate; //czêstotliwoœæ
	unsigned long byteRate; //czêstotliwoœæ bajtów, czêstotliwoœæ * liczba kana³ów * rozdzielczoœæ / 8
	unsigned short blockAlign; //rozmiar próbki, liczba kana³ów * rozdzielczoœæ / 8
	unsigned short bitsPerSample; //rozdzielczoœæ w bitach
};

struct wave_chunk
{
	char ID[4]; //ID danych, pocz¹tek czêœci z danymi
	unsigned long size; //rozmiar bloku danych
};

class WaveReader
{
private:
	wave_header header; //blok RIFF
	FILE *fin; //plik wav
	ofstream output; //plik wynikowy
	wave_chunk chunk; //blok data
	int sampleSize; //rozmiar sampla
	int sampleCount; //iloœæ sampli w pliku
	int N; //iloœæ sampli w kanale
	short int *Lchannel; //lewy kana³
	short int *Rchannel; //prawy kana³
	_int16 XavgSqr; //œrednia kwadratowa dla danych podstawowych
	_int16 EavgSqr; //œrednia kwadratowa dla danych ró¿nicowych
	_int16 Xentropy; //H(s) dla danych podstawowych
	_int16 Eentropy; //H(s) dla danych ró¿nicowych
	tm time;
public:
	WaveReader(const char *FileName)
	{
		
		fin = fopen(FileName, "rb");
		fread(&header, sizeof(header), 1, fin);
		while (true)
		{
			fread(&chunk, sizeof(chunk), 1, fin);
			if (*(unsigned int *)&chunk.ID == 0x61746164)
			{
				break;
			}
		}
		sampleSize = header.bitsPerSample / 8;
		sampleCount = chunk.size * 8 / header.bitsPerSample;
		N = sampleCount / 2;
		Lchannel = new short int[N];
		Rchannel = new short int[N];
		DisplayInformation(FileName);
		string FileNameOut;
		FileNameOut = "Result_" + to_string(time.tm_mday) + to_string(time.tm_mon) + to_string(time.tm_year) + to_string(time.tm_hour);
		output.open(FileNameOut + ".txt");
		output << "Nazwa\tŒrednia kwadratowa X\tŒrednia kwadratowa E\tEntropia X\tEntropia E\n";
		ReadData();
		FinalValues();
		WriteToFile(FileName);

	}
	 
	~WaveReader()
	{
		fclose(fin);
		output.close();
	}

	void DisplayInformation(char const *FileName) //wyœwietlanie danych
	{
		cout << "************WAV file properties************" << endl;
		cout << "File name: " << FileName << endl;
		cout << "File type: " << header.chunkID << endl;
		cout << "File format: " << header.format << endl;
		cout << "File size: " << header.chunkSize << endl;
		cout << "Format name: " << header.subchunk1ID << endl;
		cout << "Format length: " << header.subchunk1Size << endl;
		cout << "Conversion type: " << header.audioFormat << endl;
		cout << "Number of channels: " << header.numChannels << endl;
		cout << "Sample rate: " << header.sampleRate << endl;
		cout << "Bits per sample: " << header.bitsPerSample << endl;
	}

	void ReadData() //czytanie wartoœci próbek do wektorów kana³ów
	{
		int n = 0, k = 0;
		for (int i = 0; i < sampleCount; i++)
		{
			if (i % 2 == 0) //parzyste próbki to kana³ lewy, próbki s¹ w pliku u³o¿one na przemian lewy kana³ prawy kana³ itd.
			{
				fread(&Lchannel[n], sampleSize, 1, fin);
				n++;
			}
			else
			{
				fread(&Rchannel[k], sampleSize, 1, fin);
				k++;
			}
		}
	}

	short int *Differental(short int *channel) //ró¿nicowanie kana³u
	{
		short int *differentalChannel = new short int[N];
		differentalChannel[0] = channel[0];
		for (int i = 1; i < N ; i++)
		{
			differentalChannel[i] = channel[i] - channel[i - 1];
		}
		return differentalChannel;
	}

	_int16 Calculations(short int *channel, int option)
	{
		_int16 result = 0;
		_int16 avg = 0;
		switch (option)
		{
		case 1: //œrednia kwadratowa dla danych podstawowych
			for (int i = 0; i < N - 1; i++)
			{
				avg = avg + (_int16)pow(channel[i],2);
			}
			result = (avg / N);
			break;

		default:
			cout << "Niepoprawna opcja" << endl;
			result = -1;
			break;
		}

		return result;
	}

	void FinalValues() //liczenie ostatecznych wartoœci
	{
		XavgSqr = 1/2 * (Calculations(Lchannel, 1) + Calculations(Rchannel, 1));
		EavgSqr = 1 / 2 * (Calculations(Differental(Lchannel), 1) + Calculations(Differental(Rchannel), 1));
	}

	void WriteToFile(const char *FileName) // zapisywanie wyników do txt
	{
		output << FileName << "\t" << XavgSqr << "\t" << EavgSqr <</* "\t" << Xentropy << "\t" << Eentropy << */endl;
	}

};

void main()
{
	WaveReader wave("gun.wav");
	getchar();
}