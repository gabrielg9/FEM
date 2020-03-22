#include<iostream>
#include<ctime>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<iomanip>

using namespace std;

class Node
{
public:
	double x, y, t;
	bool bc = false;

	Node(double x, double y, double t);
	~Node();
};


class Element
{
	int id; 
	int *node_id = new int[4];//bo ka¿dy element ma 4 wierzcho³ki
public:
	double x[4];
	double y[4];
	double ksi[4];
	double eta[4];
	double N1[4], N2[4], N3[4], N4[4];
	double xp[4], yp[4];
	double dN_dksi[4][4];
	double dN_deta[4][4];
	double Jacobian_2d_ksi_eta[4][4];
	double detJ[4];
	double Jacobian_2d_x_y[4][4];
	double dN_dx[4][4];
	double dN_dy[4][4];
	double dN_dx_dN_dx1[4][4];
	double dN_dx_dN_dx2[4][4];
	double dN_dx_dN_dx3[4][4];
	double dN_dx_dN_dx4[4][4];
	double dN_dy_dN_dy1[4][4];
	double dN_dy_dN_dy2[4][4];
	double dN_dy_dN_dy3[4][4];
	double dN_dy_dN_dy4[4][4];
	double dN_dx_dN_dx1_detJ[4][4];
	double dN_dx_dN_dx2_detJ[4][4];
	double dN_dx_dN_dx3_detJ[4][4];
	double dN_dx_dN_dx4_detJ[4][4];
	double dN_dy_dN_dy1_detJ[4][4];
	double dN_dy_dN_dy2_detJ[4][4];
	double dN_dy_dN_dy3_detJ[4][4];
	double dN_dy_dN_dy4_detJ[4][4];
	double k1[4][4];
	double k2[4][4];
	double k3[4][4];
	double k4[4][4];
	double **macierzH;
	double N_ksi_eta[4][4];
	double tab1C[4][4];
	double tab2C[4][4];
	double tab3C[4][4];
	double tab4C[4][4];
	double **macierzC;
	double L[4];//dlugosc boku
	double DetJ, DetJ_1, DetJ_2, DetJ_3, DetJ_4;
	double pow1[4][4], pow2[4][4], pow3[4][4], pow4[4][4];
	double N_WB[2][4], N_WB_1[2][4], N_WB_2[2][4], N_WB_3[2][4], N_WB_4[2][4];
	double ksi2[2], eta2[2];
	double tab_bok[4];
	double **matrix_HBC;
	double **matrixH_final;
	double *wektorP;
	double N_pow1[2][4], N_pow2[2][4], N_pow3[2][4], N_pow4[2][4];
	



	void oblicz_x_y(Node**);
	void oblicz_ksi_eta();
	void oblicz_f_ksztaltu();
	void oblicz_pkt_calkowania();
	void pochodne_ksi_eta();
	void Jac_2d_ksi_eta();
	void wyznacznik_Jacobianu();
	void pochodna_Jacobianu();
	void dN_dxx();
	void dN_dyy();
	void skladowe();
	void skladowe_wyznacznik();
	void skladowe_wyznacznik_konduktywnosc(double);
	double **matrixH();
	void skladowe_do_mac_C_ksi_eta();
	void pkt_calkowania_1(double, double);
	void pkt_calkowania_2(double, double);
	void pkt_calkowania_3(double, double);
	void pkt_calkowania_4(double, double);
	double **matrixC();
	void dlugosc_boku();
	void warunki_brzegowe(double);
	void wyswietlenie_powierzchni_w_brzegowych();
	double **macierz_H_HBC();
	double **matrix_final();
	double *vectorP(double, double);

	Element(int id, int *node_id);
	Element();
	~Element();
	int getID();
	int *getNodeID();
};


class Grid
{
	vector<Node*> gridNode;
	Element** gridElement;
	double **globalH;
	double **globalC;
	double **globalHBC;
	double **globalHHBC;
	double *globalP;
	double *globalP_final;
public:
	Grid(double, double, double, double, double);
	~Grid();
	void showGrid();
	vector<Node*> getNodes();
	Element** getElements();
	void Global(Element**, double, double**, double*, double, double, double);
	void showGlobalMatrix(double);
	void nowaMetoda(double*, double**, double);
};

const double eps = 1e-12;
bool gauss(int n, double ** AB, double * X)
{
	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników
	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych
	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}

int main()
{
	fstream plik;
	plik.open("plik.txt");
	double initial_temperature;//temp pocz¹tkowa 100
	double simulation_time;//500
	double step_time;//50
	double t_alfa;//1200
	double alfa;//300
	double spec_heat;//700
	double H;// wysokoœæ
	double L;//szerokoœæ
	double nH;//l. wêz³ów po wys
	double nL;//l. wêz³ów po szer
	double K;//konduktywnosc 25
	double ro;//gestosc 7800
	if (plik.good() == true)
	{
		plik >> initial_temperature >> simulation_time >> step_time >> t_alfa >> alfa >> H >> L >> nH >> nL >> spec_heat >> K >> ro ;
	}
	double size = (nH - 1) * (nL - 1);
	double size2 = nH * nL;
	
	Grid siatka(H, L, nH, nL, initial_temperature);
	//siatka.showGrid();
	plik.close();

	Element** elements = new Element*[size];
	elements = siatka.getElements();
	vector<Node*> nodes = siatka.getNodes();
	Node **tempNodes = new Node*[4];
	
	for (int i = 0; i < size; i++) {
		int *nodesId = elements[i]->getNodeID();
		for (int i = 0; i < 4; i++)
		{
			tempNodes[i] = nodes[nodesId[i] - 1];
		}
		elements[i]->oblicz_x_y(tempNodes);
		elements[i]->oblicz_ksi_eta();
		elements[i]->oblicz_f_ksztaltu();
		elements[i]->oblicz_pkt_calkowania();
		elements[i]->pochodne_ksi_eta();
		elements[i]->Jac_2d_ksi_eta();
		elements[i]->wyznacznik_Jacobianu();
		elements[i]->pochodna_Jacobianu();
		elements[i]->dN_dxx();
		elements[i]->dN_dyy();
		elements[i]->skladowe();
		elements[i]->skladowe_wyznacznik();
		elements[i]->skladowe_wyznacznik_konduktywnosc(K);
		elements[i]->matrixH();
		//cout << endl << endl << "Obliczanie macierzy C" << endl << endl;
		elements[i]->skladowe_do_mac_C_ksi_eta();
		elements[i]->pkt_calkowania_1(spec_heat, ro);
		elements[i]->pkt_calkowania_2(spec_heat, ro);
		elements[i]->pkt_calkowania_3(spec_heat, ro);
		elements[i]->pkt_calkowania_4(spec_heat, ro);
		elements[i]->matrixC();
		//cout << endl << endl << "Warunki brzegowe" << endl;
		elements[i]->dlugosc_boku();

		elements[i]->DetJ_1 = elements[i]->L[0] / 2.0;
		elements[i]->DetJ_2 = elements[i]->L[1] / 2.0;
		elements[i]->DetJ_3 = elements[i]->L[2] / 2.0;
		elements[i]->DetJ_4 = elements[i]->L[3] / 2.0;

		if(tempNodes[0]->bc && tempNodes[1]->bc)
			elements[i]->tab_bok[0] = 1;
		else
			elements[i]->tab_bok[0] = 0;
		if (tempNodes[1]->bc && tempNodes[2]->bc)
			elements[i]->tab_bok[1] = 1;
		else
			elements[i]->tab_bok[1] = 0;
		if (tempNodes[2]->bc && tempNodes[3]->bc)
			elements[i]->tab_bok[2] = 1;
		else
			elements[i]->tab_bok[2] = 0;
		if (tempNodes[3]->bc && tempNodes[0]->bc)
			elements[i]->tab_bok[3] = 1;
		else
			elements[i]->tab_bok[3] = 0;
		
		elements[i]->warunki_brzegowe(alfa);
		elements[i]->wyswietlenie_powierzchni_w_brzegowych();
		elements[i]->macierz_H_HBC();
		elements[i]->matrix_final();
		elements[i]->vectorP(t_alfa, alfa);
	}

	double **AB, *X;
	int      n, i;
	cout << setprecision(4) << fixed;

	// odczytujemy liczbê niewiadomych nodow
	n = size2;

	// tworzymy macierze AB i X
	AB = new double *[n];
	X = new double[n];
	for (i = 0; i < n; i++) AB[i] = new double[n + 1];

	double t_zero[16];
	for (int i = 0; i < 16; i++)
		t_zero[i] = initial_temperature;
	
	int ilosc_iteracji = simulation_time / step_time;
	siatka.Global(elements, size, AB, t_zero, step_time, t_alfa, alfa);
	for (int k = 0; k < ilosc_iteracji; k++)
	{
		siatka.nowaMetoda(t_zero, AB, step_time);
		siatka.showGlobalMatrix(size2);

		/*for (int i = 0; i < size2; i++)
		{
			for (int j = 0; j < size2 + 1; j++)
				cout << AB[i][j] << "  ";
			cout << endl;
		}
		cout << endl << endl;*/

		cout << endl << "X - iteracja nr " << k+1 << endl;
		if (gauss(n, AB, X))
		{
			for (i = 0; i < n; i++)
				cout << "x" << i + 1 << " = " << setw(9) << X[i] << endl;
		}
		else
			cout << "DZIELNIK ZERO\n";
		for (int r = 0; r < 16; r++)
			t_zero[r] = X[r];
	}
	
	for (i = 0; i < n; i++) delete[] AB[i];
	delete[] AB;
	delete[] X;

	system("pause");
	return 0;
}

Node::Node(double x, double y, double t)
{
	this->x = x;
	this->y = y;
	this->t = t;
}
Node::~Node()
{

}
Element::Element(int id, int *node_id)
{
	this->id = id;
	this->node_id = node_id;

	macierzC = new double*[4];
	for (int i = 0; i < 4; i++)
		macierzC[i] = new double[4];

	macierzH = new double*[4];
	for (int i = 0; i < 4; i++)
		macierzH[i] = new double[4];

	matrix_HBC = new double*[4];
	for (int i = 0; i < 4; i++)
		matrix_HBC[i] = new double[4];

	matrixH_final = new double*[4];
	for (int i = 0; i < 4; i++)
		matrixH_final[i] = new double[4];

	wektorP = new double[4];
}
Element::~Element()
{

}
int Element::getID()
{
	return id;
}
int *Element::getNodeID()
{
	return node_id;
}
Grid::Grid(double H, double L, double nH, double nL, double to)
{
	double nodeHeight = H / (nH - 1);
	double nodeWidth = L / (nL - 1);
	for (int i = 0; i < nL; i++)
	{
		for (int j = 0; j < nH; j++)
		{
			Node *node = new Node(i*nodeWidth, j*nodeHeight, to);
			if ((i == 0 || i == nL - 1) || (j == nH - 1 || j == 0))
				node->bc = true;
			gridNode.push_back(node);
		}
	}
	gridElement = new Element *[9];
	int id = 0;
	for (int i = 0; i < nL - 1; i++)
	{
		for (int j = 0; j < nH - 1; j++)
		{
			int *nodeId = new int[4];
			nodeId[0] = (j + 1) + nH * i;
			nodeId[1] = (j + 1) + nH + nH * i;
			nodeId[2] = (j + 2) + nH + nH * i;
			nodeId[3] = (j + 2) + nH * i;

			gridElement[id++] = new Element(id+1, nodeId);

		}
	}
	globalH = new double *[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		globalH[i] = new double[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		for (int j = 0; j < gridNode.size(); j++)
			globalH[i][j] = 0;

	globalC = new double *[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		globalC[i] = new double[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		for (int j = 0; j < gridNode.size(); j++)
			globalC[i][j] = 0;

	globalHBC = new double *[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		globalHBC[i] = new double[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		for (int j = 0; j < gridNode.size(); j++)
			globalHBC[i][j] = 0;

	globalHHBC = new double *[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		globalHHBC[i] = new double[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		for (int j = 0; j < gridNode.size(); j++)
			globalHHBC[i][j] = 0;

	globalP = new double[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		globalP[i] = 0;

	globalP_final = new double[gridNode.size()];
	for (int i = 0; i < gridNode.size(); i++)
		globalP_final[i] = 0;
}

Grid::~Grid()
{

}

void Grid::showGrid()
{
	for (int i = 0; i < 9; i++)
	{
		cout << gridElement[i]->getID() << " = ";
		int *nodeID = gridElement[i]->getNodeID();
		cout << nodeID[0] << " " << nodeID[1] << " " << nodeID[2] << " " << nodeID[3] << endl;
	}
}

vector<Node*> Grid::getNodes()
{
	return gridNode;
}

Element** Grid::getElements()
{
	return gridElement;
}

void Grid::Global(Element** elements, double size, double** AB, double* t_zero, double step_time, double t_alfa, double alfa)
{
	for (int k = 0; k < size; k++)
	{
		globalP_final[k] = 0;
		for (int r = 0; r < size; r++)
		{
			globalH[k][r] = 0;
			globalC[k][r] = 0;
			globalHBC[k][r] = 0;
			globalHHBC[k][r] = 0;
		}
	}
	for (int i = 0; i < size; i++)
	{
		double **localH = elements[i]->matrixH();
		double **localC = elements[i]->matrixC();
		double **localHBC = elements[i]->macierz_H_HBC();
		double *localP = elements[i]->vectorP(t_alfa, alfa);
		int *nodeId = elements[i]->getNodeID();
		
		for (int j = 0; j < 4; j++)
		{
			globalP[nodeId[j] - 1] += localP[j];
			for (int k = 0; k < 4; k++)
			{
				globalH[nodeId[j] - 1][nodeId[k] - 1] += localH[j][k];
				globalC[nodeId[j] - 1][nodeId[k] - 1] += localC[j][k];
				globalHBC[nodeId[j] - 1][nodeId[k] - 1] += localHBC[j][k];
			}
		}
	}
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			globalHHBC[i][j] += globalH[i][j] + globalHBC[i][j] + (globalC[i][j] / step_time);
		}
	for (int i = 0; i < 16; i++)
		globalP_final[i] += globalP[i];

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 17; j++)
		{
			if (j == 16)
				AB[i][j] = globalP_final[i];
			else
				AB[i][j] = globalHHBC[i][j];
		}
	}		
}

void Grid::showGlobalMatrix(double size)
{
	/*cout << endl << "Global [C]" << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			cout << globalC[i][j] << " ";
		cout << endl;
	}*/

	/*cout << endl << endl << "Global [H]" << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			cout << globalH[i][j] << " ";
		cout << endl;
	}*/

	/*cout << endl << endl << "Global [H] + [C]/dT" << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			cout << globalHHBC[i][j] << " ";
		cout << endl;
	}*/

	/*cout << endl << endl << "Vector P Global" << endl;
	for (int i = 0; i < size; i++)
		cout << globalP_final[i] << "  ";*/
}

void Grid::nowaMetoda(double *t_zero, double **AB, double step_time)
{
	for (int i = 0; i < 16; i++)
		globalP_final[i] = 0;
	for (int i = 0; i < 16; i++)
		for (int j = 0; j < 16; j++)
		{
			globalP_final[i] += (globalC[i][j] / step_time)*t_zero[j];
		}
	for (int i = 0; i < 16; i++)
		globalP_final[i] += globalP[i];

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 17; j++)
		{
			if (j == 16)
				AB[i][j] = globalP_final[i];
			else
				AB[i][j] = globalHHBC[i][j];
		}
	}
}




Element::Element()
{
	
}

void Element::oblicz_x_y(Node **node)
{
	x[0] = node[0]->x; x[1] = node[1]->x; x[2] = node[2]->x; x[3] = node[3]->x;
	y[0] = node[0]->y; y[1] = node[1]->y; y[2] = node[2]->y; y[3] = node[3]->y;
}

void Element::oblicz_ksi_eta()
{
	ksi[0] = -(1 / sqrt(3));
	ksi[1] = -ksi[0];
	ksi[2] = ksi[1];
	ksi[3] = -ksi[2];

	eta[0] = ksi[0];
	eta[1] = eta[0];
	eta[2] = -eta[1];
	eta[3] = eta[2];
}

void Element::oblicz_f_ksztaltu()
{
	for (int i = 0; i < 4; i++)
	{
		N1[i] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		N2[i] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		N3[i] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		N4[i] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}

	/*for (int i = 0; i < 4; i++)
		cout << endl << "N1-N4" << endl << N1[i] << "  " << N2[i] << "  " << N3[i] << "  " << N4[i] << endl;*/
}

void Element::oblicz_pkt_calkowania()
{
	for (int i = 0; i < 4; i++)
	{
		xp[i] = N1[i] * x[0] + N2[i] * x[1] + N3[i] * x[2] + N4[i] * x[3];
		yp[i] = N1[i] * y[0] + N2[i] * y[1] + N3[i] * y[2] + N4[i] * y[3];
	}

	/*for (int i = 0; i < 4; i++)
		cout << endl << endl << xp[i] << "  " << yp[i] << endl;*/
}

void Element::pochodne_ksi_eta()
{
	for (int i = 0; i < 4; i++)
	{
		dN_dksi[0][i] = -0.25*(1 - eta[i]);
		dN_dksi[1][i] = 0.25*(1 - eta[i]);
		dN_dksi[2][i] = 0.25*(1 + eta[i]);
		dN_dksi[3][i] = -0.25*(1 + eta[i]);
	}

	/*cout << endl << endl << "dN dksi" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << dN_dksi[i][j] << "  ";
		cout << endl;
	}*/

	for (int i = 0; i < 4; i++)
	{
		dN_deta[0][i] = -0.25*(1 - ksi[i]);
		dN_deta[1][i] = -0.25*(1 + ksi[i]);
		dN_deta[2][i] = 0.25*(1 + ksi[i]);
		dN_deta[3][i] = 0.25*(1 - ksi[i]);
	}

	/*cout << endl << endl << "dN dEta" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			//cout << dN_deta[i][j] << "  ";
		//cout << endl;
	}*/
}

void Element::Jac_2d_ksi_eta()
{
	for (int i = 0; i < 4; i++)
	{
		Jacobian_2d_ksi_eta[0][i] = dN_dksi[0][i] * x[0] + dN_dksi[1][i] * x[1] + dN_dksi[2][i] * x[2] + dN_dksi[3][i] * x[3];
		Jacobian_2d_ksi_eta[1][i] = dN_dksi[0][i] * y[0] + dN_dksi[1][i] * y[1] + dN_dksi[2][i] * y[2] + dN_dksi[3][i] * y[3];
		Jacobian_2d_ksi_eta[2][i] = dN_deta[0][i] * x[0] + dN_deta[1][i] * x[1] + dN_deta[2][i] * x[2] + dN_deta[3][i] * x[3];
		Jacobian_2d_ksi_eta[3][i] = dN_deta[0][i] * y[0] + dN_deta[1][i] * y[1] + dN_deta[2][i] * y[2] + dN_deta[3][i] * y[3];
	}

	/*cout << endl << endl << "Jacobian" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << Jacobian_2d_ksi_eta[i][j] << "  ";
		cout << endl;
	}*/
}

void Element::wyznacznik_Jacobianu()
{
	for (int i = 0; i < 4; i++)
		detJ[i] = Jacobian_2d_ksi_eta[0][i] * Jacobian_2d_ksi_eta[3][i] - Jacobian_2d_ksi_eta[1][i] * Jacobian_2d_ksi_eta[2][i];

	/*cout << endl << endl << "Wyznacznik Jacobianu" << endl;
	for (int i = 0; i < 4; i++)
		cout << detJ[i] << "  ";*/
}

void Element::pochodna_Jacobianu()
{
	for (int i = 0; i < 4; i++)
	{
		Jacobian_2d_x_y[0][i] = Jacobian_2d_ksi_eta[3][i] / detJ[i];
		Jacobian_2d_x_y[1][i] = Jacobian_2d_ksi_eta[1][i] / detJ[i];
		Jacobian_2d_x_y[2][i] = Jacobian_2d_ksi_eta[2][i] / detJ[i];
		Jacobian_2d_x_y[3][i] = Jacobian_2d_ksi_eta[0][i] / detJ[i];
	}

	/*cout << endl << endl << "Pochodna Jacobianu" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << Jacobian_2d_x_y[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::dN_dxx()
{
	for (int i = 0; i < 4; i++)
	{
		dN_dx[0][i] = Jacobian_2d_x_y[0][0] * dN_dksi[i][0] + Jacobian_2d_x_y[1][0] * dN_deta[i][0];
		dN_dx[1][i] = Jacobian_2d_x_y[0][1] * dN_dksi[i][1] + Jacobian_2d_x_y[1][1] * dN_deta[i][1];
		dN_dx[2][i] = Jacobian_2d_x_y[0][2] * dN_dksi[i][2] + Jacobian_2d_x_y[1][2] * dN_deta[i][2];
		dN_dx[3][i] = Jacobian_2d_x_y[0][3] * dN_dksi[i][3] + Jacobian_2d_x_y[1][3] * dN_deta[i][3];
	}

	/*cout << endl << endl << "dN dx" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << dN_dx[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::dN_dyy()
{
	for (int i = 0; i < 4; i++)
	{
		dN_dy[0][i] = Jacobian_2d_x_y[2][0] * dN_dksi[i][0] + Jacobian_2d_x_y[3][0] * dN_deta[i][0];
		dN_dy[1][i] = Jacobian_2d_x_y[2][1] * dN_dksi[i][1] + Jacobian_2d_x_y[3][1] * dN_deta[i][1];
		dN_dy[2][i] = Jacobian_2d_x_y[2][2] * dN_dksi[i][2] + Jacobian_2d_x_y[3][2] * dN_deta[i][2];
		dN_dy[3][i] = Jacobian_2d_x_y[2][3] * dN_dksi[i][3] + Jacobian_2d_x_y[3][3] * dN_deta[i][3];
	}

	/*cout << endl << endl << "dN dy" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << dN_dy[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::skladowe()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dN_dx_dN_dx1[j][i] = dN_dx[0][j] * dN_dx[0][i];
			dN_dx_dN_dx2[j][i] = dN_dx[1][j] * dN_dx[1][i];
			dN_dx_dN_dx3[j][i] = dN_dx[2][j] * dN_dx[2][i];
			dN_dx_dN_dx4[j][i] = dN_dx[3][j] * dN_dx[3][i];
			dN_dy_dN_dy1[j][i] = dN_dy[0][j] * dN_dy[0][i];
			dN_dy_dN_dy2[j][i] = dN_dy[1][j] * dN_dy[1][i];
			dN_dy_dN_dy3[j][i] = dN_dy[2][j] * dN_dy[2][i];
			dN_dy_dN_dy4[j][i] = dN_dy[3][j] * dN_dy[3][i];
		}
	}

	/*cout << endl << endl << "dN dx dN dx lub dN dy dN dy" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << dN_dy_dN_dy4[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::skladowe_wyznacznik()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dN_dx_dN_dx1_detJ[i][j] = dN_dx_dN_dx1[i][j] * detJ[0];
			dN_dx_dN_dx2_detJ[i][j] = dN_dx_dN_dx2[i][j] * detJ[1];
			dN_dx_dN_dx3_detJ[i][j] = dN_dx_dN_dx3[i][j] * detJ[2];
			dN_dx_dN_dx4_detJ[i][j] = dN_dx_dN_dx4[i][j] * detJ[3];
			dN_dy_dN_dy1_detJ[i][j] = dN_dy_dN_dy1[i][j] * detJ[0];
			dN_dy_dN_dy2_detJ[i][j] = dN_dy_dN_dy2[i][j] * detJ[1];
			dN_dy_dN_dy3_detJ[i][j] = dN_dy_dN_dy3[i][j] * detJ[2];
			dN_dy_dN_dy4_detJ[i][j] = dN_dy_dN_dy4[i][j] * detJ[3];
		}
	}

	/*cout << endl << endl << "Pomnozone razy wyznacznik" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << dN_dy_dN_dy4_detJ[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::skladowe_wyznacznik_konduktywnosc(double k)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			k1[i][j] = k * (dN_dx_dN_dx1_detJ[i][j] + dN_dy_dN_dy1_detJ[i][j]);
			k2[i][j] = k * (dN_dx_dN_dx2_detJ[i][j] + dN_dy_dN_dy2_detJ[i][j]);
			k3[i][j] = k * (dN_dx_dN_dx3_detJ[i][j] + dN_dy_dN_dy3_detJ[i][j]);
			k4[i][j] = k * (dN_dx_dN_dx4_detJ[i][j] + dN_dy_dN_dy4_detJ[i][j]);
		}
	}

	/*cout << endl << endl << "Pomnozone razy konduktywnosc" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << k4[i][j] << "   ";
		cout << endl;
	}*/
}

double ** Element::matrixH()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			macierzH[i][j] = k1[i][j] + k2[i][j] + k3[i][j] + k4[i][j];
	}

	double **tab3 = new double*[4];
	for (int i = 0; i < 4; i++)
		tab3[i] = macierzH[i];

	/*cout << endl << endl << "Macierz H " << endl;;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << macierzH[i][j] << "   ";
		cout << endl;
	}*/
	return tab3;
}

void Element::skladowe_do_mac_C_ksi_eta()
{
	for (int i = 0; i < 4; i++)
	{
		N_ksi_eta[i][0] = (1 - ksi[i])*(1 - eta[i]) / 4;
		N_ksi_eta[i][1] = (1 + ksi[i])*(1 - eta[i]) / 4;
		N_ksi_eta[i][2] = (1 + ksi[i])*(1 + eta[i]) / 4;
		N_ksi_eta[i][3] = (1 - ksi[i])*(1 + eta[i]) / 4;
	}
	/*cout << "Skladowe N1...N4 wzgledem ksi i eta" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << N_ksi_eta[i][j] << "    ";
		cout << endl;
	}*/
}

void Element::pkt_calkowania_1(double c, double ro)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j<4; j++)
			tab1C[i][j] = N_ksi_eta[0][j] * N_ksi_eta[i][0] * detJ[0] * c*ro;
	}

	/*cout << endl << endl << "1 Punkt calkowania" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << tab1C[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::pkt_calkowania_2(double c, double ro)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j<4; j++)
			tab2C[i][j] = N_ksi_eta[1][j] * N_ksi_eta[i][1] * detJ[1] * c*ro;
	}

	/*cout << endl << endl << "2 Punkt calkowania" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << tab2C[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::pkt_calkowania_3(double c, double ro)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j<4; j++)
			tab3C[i][j] = N_ksi_eta[2][j] * N_ksi_eta[i][2] * detJ[2] * c*ro;
	}

	/*cout << endl << endl << "3 Punkt calkowania" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << tab3C[i][j] << "   ";
		cout << endl;
	}*/
}

void Element::pkt_calkowania_4(double c, double ro)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j<4; j++)
			tab4C[i][j] = N_ksi_eta[3][j] * N_ksi_eta[i][3] * detJ[3] * c*ro;
	}

	/*cout << endl << endl << "4 Punkt calkowania" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << tab4C[i][j] << "   ";
		cout << endl;
	}*/
}

double ** Element::matrixC()
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			macierzC[i][j] = (tab1C[i][j] + tab2C[i][j] + tab3C[i][j] + tab4C[i][j]);
	double **tab1 = new double*[4];
	for (int i = 0; i < 4; i++)
		tab1[i] = macierzC[i];

	/*cout << endl << endl << "Macierz C" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << macierzC[i][j] << "    ";
		cout << endl;
	}*/
	return tab1;
}

void Element::dlugosc_boku()
{
	for (int i = 0; i<4; i++)
		L[i] = sqrt(pow((x[(i + 1) % 4] - x[i]), 2) + pow((y[(i + 1) % 4] - y[i]), 2));

	/*cout << endl << endl << "Dlugosci boku" << endl;
	for (int i = 0; i < 4; i++)
		cout << L[i] << "    ";
	cout << endl;*/
}

void Element::warunki_brzegowe(double alfa)
{
	for (int i = 0; i < 4; i++)
	{
		DetJ = L[i] / 2.0;
		if (i == 0)
		{
			ksi2[0] = -1 / sqrt(3);
			ksi2[1] = -ksi2[0];
			eta2[0] = -1;
			eta2[1] = -1;
		}
		else if (i == 1)
		{
			ksi2[0] = 1;
			ksi2[1] = 1;
			eta2[0] = -1 / sqrt(3);
			eta2[1] = -eta2[0];
		}
		else if (i == 2)
		{
			ksi2[0] = 1 / sqrt(3);
			ksi2[1] = -ksi2[0];
			eta2[0] = 1;
			eta2[1] = 1;
		}
		else
		{
			ksi2[0] = -1;
			ksi2[1] = -1;
			eta2[0] = 1 / sqrt(3);
			eta2[1] = -eta2[0];
		}

		for (int j = 0; j < 2; j++)
		{
			N_WB[j][0] = (1 - ksi2[j])*(1 - eta2[j]) / 4;
			N_WB[j][1] = (1 + ksi2[j])*(1 - eta2[j]) / 4;
			N_WB[j][2] = (1 + ksi2[j])*(1 + eta2[j]) / 4;
			N_WB[j][3] = (1 - ksi2[j])*(1 + eta2[j]) / 4;
			if (i == 0)
			{
				N_WB_1[j][0] = N_WB[j][0];
				N_WB_1[j][1] = N_WB[j][1];
				N_WB_1[j][2] = N_WB[j][2];
				N_WB_1[j][3] = N_WB[j][3];
			}
			else if (i == 1)
			{
				N_WB_2[j][0] = N_WB[j][0];
				N_WB_2[j][1] = N_WB[j][1];
				N_WB_2[j][2] = N_WB[j][2];
				N_WB_2[j][3] = N_WB[j][3];
			}
			else if (i == 2)
			{
				N_WB_3[j][0] = N_WB[j][0];
				N_WB_3[j][1] = N_WB[j][1];
				N_WB_3[j][2] = N_WB[j][2];
				N_WB_3[j][3] = N_WB[j][3];
			}
			else
			{
				N_WB_4[j][0] = N_WB[j][0];
				N_WB_4[j][1] = N_WB[j][1];
				N_WB_4[j][2] = N_WB[j][2];
				N_WB_4[j][3] = N_WB[j][3];
			}
		}
		double pc1[4][4], pc2[4][4];
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				pc1[i][j] = N_WB[0][j] * N_WB[0][i] * alfa;
				pc2[i][j] = N_WB[1][j] * N_WB[1][i] * alfa;
			}
		}
		if (i == 0)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
					pow1[i][j] = (pc1[i][j] + pc2[i][j])* DetJ;
			}
		}
		else if (i == 1)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
					pow2[i][j] = (pc1[i][j] + pc2[i][j])* DetJ;
			}
		}
		else if (i == 2)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
					pow3[i][j] = (pc1[i][j] + pc2[i][j])* DetJ;
			}
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
					pow4[i][j] = (pc1[i][j] + pc2[i][j])* DetJ;
			}
		}

	}
}

void Element::wyswietlenie_powierzchni_w_brzegowych()
{
	//cout << endl << endl << " 4 powierzchnie" << endl;
	/*for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << pow1[i][j] << "    ";
		cout << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << pow2[i][j] << "    ";
		cout << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << pow3[i][j] << "    ";
		cout << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << pow4[i][j] << "    ";
		cout << endl;
	}*/
}

double ** Element::macierz_H_HBC()
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrix_HBC[i][j] = tab_bok[0] * pow1[i][j] + tab_bok[1] * pow2[i][j] + tab_bok[2] * pow3[i][j] + tab_bok[3] * pow4[i][j];

	/*cout << endl << endl << "Matrix HBC" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << matrix_HBC[i][j] << " ";
		cout << endl;
	}*/
	double ** tab = new double*[4];
	for (int i = 0; i < 4; i++)
		tab[i] = matrix_HBC[i];
	return tab;
}

double ** Element::matrix_final()
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrixH_final[i][j] = (macierzH[i][j] + matrix_HBC[i][j]);

	/*cout << endl << "Matrix H final" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << matrixH_final[i][j] << " ";
		cout << endl;
	}*/
	
	double ** tab = new double* [4];
	for (int i = 0; i < 4; i++)
		tab[i] = matrixH_final[i];
	return tab;
}

double * Element::vectorP(double t_alfa, double alfa)
{
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			N_pow1[i][j] = alfa * N_WB_1[i][j] * t_alfa * DetJ_1 * tab_bok[0];
			N_pow2[i][j] = alfa * N_WB_2[i][j] * t_alfa * DetJ_2 * tab_bok[1];
			N_pow3[i][j] = alfa * N_WB_3[i][j] * t_alfa * DetJ_3 * tab_bok[2];
			N_pow4[i][j] = alfa * N_WB_4[i][j] * t_alfa * DetJ_4 * tab_bok[3];
		}
	}

	for (int i = 0; i < 4; i++)
	{
		wektorP[i] = ((N_pow1[0][i] + N_pow1[1][i]) + (N_pow2[0][i] + N_pow2[1][i]) + (N_pow3[0][i] + N_pow3[1][i]) + (N_pow4[0][i] + N_pow4[1][i]));
	}

	/*cout << endl << " Vector P" << endl;
	for (int i = 0; i < 4; i++)
		cout << wektorP[i] << "  ";
	cout << endl << endl;*/
	double *tab2 = new double[4];
	for (int i = 0; i < 4; i++)
		tab2[i] = wektorP[i];
	return tab2;
}

