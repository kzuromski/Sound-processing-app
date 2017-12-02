#include "Wave.h"

void start(int r) {
	Wave wave1("ATrain.wav", r);
	Wave wave2("BeautySlept.wav", r);
	Wave wave3("death2.wav", r);
	Wave wave4("experiencia.wav", r);
	Wave wave5("chanchan.wav", r);
	Wave wave6("female_speech.wav", r);
	Wave wave7("FloorEssence.wav", r);
	Wave wave8("ItCouldBeSweet.wav", r);
	Wave wave9("Layla.wav", r);
	Wave wave10("LifeShatters.wav", r);
	Wave wave11("macabre.wav", r);
	Wave wave12("male_speech.wav", r);
	Wave wave13("SinceAlways.wav", r);
	Wave wave14("thear1.wav", r);
	Wave wave15("TomsDiner.wav", r);
	Wave wave16("velvet.wav", r);

	double tab[16];
	double average = 0;

	tab[0] = wave1.getAverage();
	tab[1] = wave2.getAverage();
	tab[2] = wave3.getAverage();
	tab[3] = wave4.getAverage();
	tab[4] = wave5.getAverage();
	tab[5] = wave6.getAverage();
	tab[6] = wave7.getAverage();
	tab[7] = wave8.getAverage();
	tab[8] = wave9.getAverage();
	tab[9] = wave10.getAverage();
	tab[10] = wave11.getAverage();
	tab[11] = wave12.getAverage();
	tab[12] = wave13.getAverage();
	tab[13] = wave14.getAverage();
	tab[14] = wave15.getAverage();
	tab[15] = wave16.getAverage();

	for (int i = 0; i < 16; i++)
		average += tab[i];

	average /= 16;
	ofstream averageFile("average.txt", ios::app);
	averageFile << r << " " << average << endl;
}

void main() {

	//for (int r = 10; r <= 600; r+=10)
	start(40);
	/*start(3);
	start(4);
	start(5);
	start(6);
	start(8);
	start(10);
	start(12);
	start(15);
	start(20);
	start(24);
	start(30);
	start(40);
	start(60);
	start(120);*/
	
	system("pause");
}