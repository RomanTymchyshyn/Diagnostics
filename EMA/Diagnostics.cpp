#include "Diagnostics.h"


///first stable version of diagnostics
///TODO: rewrite this file in better manner

class Patiente
{
	public:
		int n;
		vector<vector<double>> arr;
		Patiente(string nameOfFile = "")
		{
			ifstream fi(nameOfFile);
			fi >> n;
			arr = vector<vector<double>>(n, vector<double>(16, 0));
			for (int i = 0; i < n; ++i)
			{
				double t = 0;
				fi >> t;
				for (int j = 1; j <= 15; ++j)
				{
					double temp = 0;
					fi >> temp;
					arr[i][j] = temp;
				}
			}
		}
};

vector<vector<double>> H11(25, vector<double>(8));
vector<vector<double>> H12(25, vector<double>(8));
vector<vector<double>> H21(25, vector<double>(8));
vector<vector<double>> H22(12, vector<double>(8));

vector<vector<double>> Plist;
vector<vector<double>> Plist1;
vector<vector<double>> Dlist;
vector<vector<double>> Dlist1;

vector<Ellipsoid2D> Ellipsets;
vector<Ellipsoid2D> Ellipsets1;
vector<Ellipsoid2D> EllipseTS;
vector<Ellipsoid2D> EllipseTS1;

string GiveString(int i)
{
	char* temp = new char[5];
	string number;
	itoa(i, temp, 10);
	number += temp;
	return number;
}

double CalculateStatistic(vector<double> x, vector<double> y)
{
	sort(x.begin(), x.end());
	double count = 0;
	double n = x.size();
	double m = y.size();
	double N = n * (n - 1) / 2;
	double g = 3.0;
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
		{
			double h = 0.0;
			for (int k = 0; k < m; k++)
				if (x[i] < y[k] && y[k] <= x[j])
					h++;
			h /= m;
			double p1 = (h * m + g * g / 2 - g * sqrt(h * (1 - h) * m + g * g / 4)) / (m + g * g);
			double p2 = (h * m + g * g / 2 + g * sqrt(h * (1 - h) * m + g * g / 4)) / (m + g * g);
			double p = (j - i) / (n + 1.0);
			if (p1 < p && p <= p2)
				count++;
		}
	return count / N;
}

void CalculateAverageStatistics(vector<Patiente> patient, vector<Patiente> patient1)
{
	for (int k = 0; k < patient.size(); k++)
	{
		vector<double> p(15, 0);
		for (int i = 1; i <= 15; i++)
		{
			vector<double> x(patient[k].n);
			for (int j = 0; j < x.size(); j++)
				x[j] = patient[k].arr[j][i];
			for (int l = 0; l < patient.size(); l++)
			{
				vector<double> y(patient[l].n);
				if (l == k)
					continue;
				for (int j = 0; j < y.size(); j++)
					y[j] = patient[l].arr[j][i];
				p[i - 1] += CalculateStatistic(x, y);
			}
			p[i - 1] /= (patient.size() - 1);
		}
		Plist.push_back(p);
	}

	for (int k = 0; k < patient1.size(); k++)
	{
		vector<double> p1(15, 0);
		for (int i = 1; i <= 15; i++)
		{
			vector<double> x(patient1[k].n);
			for (int j = 0; j < x.size(); j++)
				x[j] = patient1[k].arr[j][i];
			for (int l = 0; l < patient.size(); l++)
			{
				vector<double> y(patient[l].n);
				for (int j = 0; j < y.size(); j++)
					y[j] = patient[l].arr[j][i];
				p1[i - 1] += CalculateStatistic(x, y);
			}
			p1[i - 1] /= patient.size();
		}
		Plist1.push_back(p1);
	}

	for (int k = 0; k < patient1.size(); k++)
	{
		vector<double> d(15, 0);
		for (int i = 1; i <= 15; i++)
		{
			vector<double> x(patient1[k].n);
			for (int j = 0; j < x.size(); j++)
				x[j] = patient1[k].arr[j][i];
			for (int l = 0; l < patient1.size(); l++)
			{
				if (l == k)
					continue;
				vector<double> y(patient1[l].n);
				for (int j = 0; j < y.size(); j++)
					y[j] = patient1[l].arr[j][i];
				d[i - 1] += CalculateStatistic(x, y);
			}
			d[i - 1] /= (patient1.size() - 1);
		}
		Dlist.push_back(d);
	}

	for (int k = 0; k < patient.size(); k++)
	{
		vector<double> d1(15, 0);
		for (int i = 1; i <= 15; i++)
		{
			vector<double> x(patient[k].n);
			for (int j = 0; j < x.size(); j++)
				x[j] = patient[k].arr[j][i];
			for (int l = 0; l < patient1.size(); l++)
			{
				vector<double> y(patient1[l].n);
				for (int j = 0; j < y.size(); j++)
					y[j] = patient1[l].arr[j][i];
				d1[i - 1] += CalculateStatistic(x, y);
			}
			d1[i - 1] /= patient1.size();
		}
		Dlist1.push_back(d1);
	}
}

void CalculateEllipse(IEllipsoidBuilder* builder)
{
	for (int t = 0; t < 15; t++)
		for (int s = t + 1; s < 15; s++)
		{
			vector<Point2D> ListP;
			for (int i = 0; i < Plist.size(); i++)
			{
				Point2D tmp = new Point2D(Plist[i][t], Plist[i][s]);
				ListP.push_back(tmp);
			}
			Ellipsets.push_back(*builder->Exec(ListP));
		}

	for (int t = 0; t < 15; t++)
		for (int s = t + 1; s < 15; s++)
		{
			vector<Point2D> ListP1;
			for (int i = 0; i < Plist1.size(); i++)
			{
				Point2D tmp = new Point2D(Plist1[i][t], Plist1[i][s]);
				ListP1.push_back(tmp);
			}
			Ellipsets1.push_back(*builder->Exec(ListP1));
		}

	for (int t = 0; t < 15; t++)
		for (int s = t + 1; s < 15; s++)
		{
			vector<Point2D> ListD;
			for (int i = 0; i < Dlist.size(); i++)
			{
				Point2D tmp = new Point2D(Dlist[i][t], Dlist[i][s]);
				ListD.push_back(tmp);
			}
			EllipseTS.push_back(*builder->Exec(ListD));
		}

	for (int t = 0; t < 15; t++)
		for (int s = t + 1; s < 15; s++)
		{
			vector<Point2D> ListD1;
			for (int i = 0; i < Dlist1.size(); i++)
			{
				Point2D tmp = new Point2D(Dlist1[i][t], Dlist1[i][s]);
				ListD1.push_back(tmp);
			}
			EllipseTS1.push_back(*builder->Exec(ListD1));
		}
}

vector<double> calculateH(vector<Patiente> patient, vector<Patiente> patient1, Patiente test)
{
	vector<double> p(15);
	for (int i = 1; i <= 15; i++)
	{
		vector<double> x(test.n);
		for (int j = 0; j < x.size(); j++)
			x[j] = test.arr[j][i];
		for (int l = 0; l < patient.size(); l++)
		{
			vector<double> y(patient[l].n);
			for (int j = 0; j < y.size(); j++)
				y[j] = patient[l].arr[j][i];
			p[i - 1] += CalculateStatistic(x, y);
		}
		p[i - 1] /= patient.size();
	}

	vector<double> d(15);
	for (int i = 1; i <= 15; i++)
	{
		vector<double> x(test.n);
		for (int j = 0; j < x.size(); j++)
			x[j] = test.arr[j][i];

		for (int l = 0; l < patient1.size(); l++)
		{
			vector<double> y(patient1[l].n);
			for (int j = 0; j < y.size(); j++)
				y[j] = patient1[l].arr[j][i];
			d[i - 1] += CalculateStatistic(x, y);
		}
		d[i - 1] /= patient1.size();
	}
	vector<double> C(4);
	bool A1, A2, A3, A4;
	bool Aa1, Aa2, Aa3, Aa4;

	int iter = 0;
	for (int t = 0; t < 15; t++)
		for (int s = t + 1; s < 15; s++)
		{
			Point2D tmp1(p[t], p[s]);
			A1 = Ellipsets[iter].Inside(tmp1);
			A2 = Ellipsets1[iter].Inside(tmp1);
			A3 = (Ellipsets[iter].Inside(tmp1) && !Ellipsets1[iter].Inside(tmp1));
			A4 = (Ellipsets1[iter].Inside(tmp1) && !Ellipsets[iter].Inside(tmp1));

			if (A3 == true)
				C[0]++;
			if (A4 == true)
				C[1]++;
			if (A1 == true)
				C[2]++;
			if (A2 == true)
				C[3]++;
			iter++;
		}

	int iter1 = 0;
	for (int t = 0; t < 15; t++)
		for (int s = t + 1; s < 15; s++)
		{
			Point2D tmp2(d[t], d[s]);
			Aa1 = EllipseTS[iter1].Inside(tmp2);
			Aa2 = EllipseTS1[iter1].Inside(tmp2);
			Aa3 = (EllipseTS[iter1].Inside(tmp2) && !EllipseTS1[iter1].Inside(tmp2));
			Aa4 = (EllipseTS1[iter1].Inside(tmp2) && !EllipseTS[iter1].Inside(tmp2));

			if (Aa4 == true)
				C[0]++;
			if (Aa3 == true)
				C[1]++;
			if (Aa2 == true)
				C[2]++;
			if (Aa1 == true)
				C[3]++;
			iter1++;
		}

	for (int i = 0; i < 4; i++)
		C[i] /= 210;
	return C;
}

void Clear()
{
	Plist.clear();
	Plist1.clear();
	Dlist.clear();
	Dlist1.clear();
	Ellipsets.clear();
	Ellipsets1.clear();
	EllipseTS.clear();
	EllipseTS1.clear();
}

void TestPatientes1(IEllipsoidBuilder* builder)
{
	ofstream fo("ResultsGenetic (25-25).txt");

	vector<Patiente> patient;
	vector<Patiente> patient1;
	for (int i = 0; i < 25; i++)
	{
		string number = GiveString(301 + i);
		patient1.push_back(Patiente("Resources/D2D3/D3/D" + number + ".POK"));
	}
	for (int k = 0; k < 25; k++)
	{
		string number = GiveString(201 + k);
		Patiente test = Patiente("Resources/D2D3/D2/D" + number + ".POK");
		patient.clear();
		for (int i = 0; i < 25; i++)
		{
			if (i != k)
			{
				string number = GiveString(201 + i);
				patient.push_back(Patiente("Resources/D2D3/D2/D" + number + ".POK"));
			}
		}
		CalculateAverageStatistics(patient, patient1);
		CalculateEllipse(builder);
		vector<double> TempC = calculateH(patient, patient1, test);
		for (int i = 0; i < TempC.size(); i++)
			H11[k][i] = TempC[i];
		
		Clear();
	}


	vector<Patiente> patientd;
	vector<Patiente> patientd1;
	for (int i = 0; i < 25; i++)
	{
		string number = GiveString(201 + i);
		patientd.push_back(Patiente("Resources/D2D3/D2/D" + number + ".POK"));
	}
	for (int k = 0; k < 25; k++)
	{
		string number = GiveString(301 + k);
		Patiente test = Patiente("Resources/D2D3/D3/D" + number + ".POK");
		patientd1.clear();
		for (int i = 0; i < 25; i++)
		{
			if (i != k)
			{
				string number = GiveString(301 + i);
				patientd1.push_back(Patiente("Resources/D2D3/D3/D" + number + ".POK"));
			}
		}
		CalculateAverageStatistics(patientd, patientd1);
		CalculateEllipse(builder);
		vector<double> TempC = calculateH(patientd, patientd1, test);
		for (int i = 0; i < TempC.size(); i++)
			H12[k][i] = TempC[i];

		Clear();
	}


	double v11 = 0;
	double v12 = 0;
	double v21 = 0;
	double v22 = 0;
	for (int k = 0; k < 25; k++)
	{
		if (H11[k][3] >= H11[k][2])
			v12++;
		else
			v11++;
		if (H12[k][3] >= H12[k][2])
			v22++;
		else
			v21++;
	}
	fo << "v11: " << v11 / 25.0 << " v12: " << v12 / 25.0 << " v22: " << v22 / 25.0 << "v21: " << v21 / 25.0 << endl;
}

void TestPatientes2(IEllipsoidBuilder* builder)
{
	ofstream fo("ResultsGenetic (25-12).txt");

	vector<Patiente> Spatient;
	vector<Patiente> Spatient1;
	for (int i = 0; i < 12; i++)
	{
		string number = GiveString(301 + i);
		Spatient1.push_back(Patiente("Resources/D2D3/D3/D" + number + ".POK"));
	}
	for (int k = 0; k < 25; k++)
	{
		string number = GiveString(201 + k);
		Patiente test = Patiente("Resources/D2D3/D2/D" + number + ".POK");
		Spatient.clear();
		for (int i = 0; i < 25; i++)
		{
			if (i != k)
			{
				string number = GiveString(201 + i);
				Spatient.push_back(Patiente("Resources/D2D3/D2/D" + number + ".POK"));
			}
		}
		CalculateAverageStatistics(Spatient, Spatient1);
		CalculateEllipse(builder);
		vector<double> TempC = calculateH(Spatient, Spatient1, test);
		for (int i = 0; i < TempC.size(); i++)
			H21[k][i] = TempC[i];
			
		Clear();
	}
	
	vector<Patiente> Spatientd;
	vector<Patiente> Spatientd1;
	for (int i = 0; i < 25; i++)
	{
		string number = GiveString(201 + i);
		Spatientd.push_back(Patiente("Resources/D2D3/D2/D" + number + ".POK"));
	}
	for (int k = 0; k < 12; k++)
	{
		string number = GiveString(301 + k);
		Patiente test = Patiente("Resources/D2D3/D3/D" + number + ".POK");
		Spatientd1.clear();
		for (int i = 0; i < 12; i++)
		{
			if (i != k)
			{
				string number = GiveString(301 + i);
				Spatientd1.push_back(Patiente("Resources/D2D3/D3/D" + number + ".POK"));
			}
		}
		CalculateAverageStatistics(Spatientd, Spatientd1);
		CalculateEllipse(builder);
		vector<double> TempC = calculateH(Spatientd, Spatientd1, test);
		for (int i = 0; i < TempC.size(); i++)
			H22[k][i] = TempC[i];
	
		Clear();
	}
	
	double v11 = 0;
	double v12 = 0;
	double v21 = 0;
	double v22 = 0;
	for (int k = 0; k < 25; k++)
	{
		if (H21[k][2] <= H21[k][3])
			v12++;
		else
			v11++;
	}
	for (int k = 0; k < 12; k++)
	{
		if (H22[k][2] <= H22[k][3])
			v22++;
		else
			v21++;
	}
	fo << "v11: " << v11 / 25.0 << " v12: " << v12 / 25.0 << " v22: " << v22 / 12.0 << "v21: " << v21 / 12.0 << endl;
}