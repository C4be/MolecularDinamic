#include <fstream>
#include <cmath>

#include "constants.h"
#include "start_cond.h"

// начальная инициализация векторов сил
inline void reset_Fxyz() { 
	for (int i = 0; i < NUMBERPARTICLES; ++i) {
		Fx[i] = .0; Fy[i] = .0; Fz[i] = .0;
	}
}

// Выделение памяти	
void memory_allocation()
{	
	coordx = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	coordy = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	coordz = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	Fx = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	Fy = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	Fz = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	reset_Fxyz();
	vx = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	vy = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
	vz = (double*)malloc(sizeof(double) * NUMBERPARTICLES);
}

// высвобождение памяти
void memory_free()
{	
	free(coordx);
	free(coordy);
	free(coordz);
	free(Fx);
	free(Fy);
	free(Fz);
	free(vx);
	free(vy);
	free(vz);
}

// вычетание второй частицы из первой
void vec_rij(int first, int second)
{
	rij[0] = coordx[first] - coordx[second];
	rij[1] = coordy[first] - coordy[second];
	rij[2] = coordz[first] - coordz[second];
}

// вычисление евклидового пространства
inline double vec_euklid(double* vec) { return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]; }

// вычисление абсолютного значения вектора
inline void vec_absolute() { rij_abs = sqrt(pow(rij[0], 2) + pow(rij[1], 2) + pow(rij[2], 2)); }

// вычисление нормализированного вектора
void vec_normalize() { 
	rij_normalize[0] = rij[0] / rij_abs;
	rij_normalize[1] = rij[1] / rij_abs;
	rij_normalize[2] = rij[2] / rij_abs;
}

// подсчет силы и потенциала взаимодействия, между двумя частицами
void calc_particle_interaction(int first, int second)
{
	vec_rij(first, second); // вычисление вектора rij
	vec_absolute(); // вычисление расстояния между частицами

	if (rij_abs > RCUT)
	{	// не учитываем взаимодействия между столь удаленными частицами
		U = .0;
		F = .0;
	}
	else
	{	// вычисление по потенциалу Леннарда-Джонса
		double t1 = SIGMA / rij_abs;
		double t2 = t1 * t1 * t1 * t1 * t1 * t1;
		t1 = t2 * t2;
		U = 4 * EPS * (t1 - t2) - UCUT; 
		F = 24 * EPS * (2 * t1 - t2) / rij_abs;
	}
}

// подсчет силы и потенциала взаимодействия, между двумя частицами
std::pair<double, double>* calc_particle_interaction(double* vec_rij)
{
	double vec_tmp_abs = sqrt(vec_euklid(vec_rij));	
	// вычисление по потенциалу Леннарда-Джонса
	double t1 = SIGMA / vec_tmp_abs;
	double t2 = t1 * t1 * t1 * t1 * t1 * t1;
	t1 = t2 * t2;
	return new std::pair<double, double>{ 4 * EPS * (t1 - t2) - UCUT, 24 * EPS * (2 * t1 - t2) / vec_tmp_abs };
}

void initial_vectors_F_virt()
{

	Epot = 0; // потенциальная энергия

	// Перерасчет сил взаимодействия
	double tx, ty, tz;
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		tx = 0; ty = 0; tz = 0;

		// перерасчет сил взаимодействия
		for (int j = 0; j < NUMBERPARTICLES; ++j)
		{
			for (int dx = -1; dx < 2; ++dx) {
				for (int dy = -1; dy < 2; ++dy) {
					for (int dz = -1; dz < 2; ++dz) {
						if (i == j && dx == 0 && dy == 0 && dz == 0) continue;
						double* vec_tmp = new double[3] {coordx[i] - (coordx[j] + dx * LX), coordy[i] - (coordy[j] + dy * LY), coordz[i] - (coordz[j] + dz * LZ)};
						double vec_tmp_abs = sqrt(vec_euklid(vec_tmp));
						std::pair<double, double>* u_and_f = calc_particle_interaction(vec_tmp); // расчет потенциала им силы взаимодействия
						double tmp_F = u_and_f->second;
						if (vec_tmp_abs > RCUT) {
							delete u_and_f;
							delete[] vec_tmp;
							continue;
						}
						tx += tmp_F * vec_tmp[0] / vec_tmp_abs;
						ty += tmp_F * vec_tmp[1] / vec_tmp_abs;
						tz += tmp_F * vec_tmp[2] / vec_tmp_abs;

						Epot += u_and_f->first; // перерасчет потенциальной энергии системы

						delete u_and_f;
						delete[] vec_tmp;
					}

				}

			}

		}
		Fx[i] = tx;
		Fy[i] = ty;
		Fz[i] = tz;
	}

	Epot /= 2; // Избавление от повторяющихся значений
}

void initial_vectors_F()
{	
	Epot = 0; // Обнуление потенциальной энергии
	double tx, ty, tz;
	for (int i = 0; i < NUMBERPARTICLES - 1; ++i) 
	{
		for (int j = i + 1; j < NUMBERPARTICLES; ++j) {
			calc_particle_interaction(i, j); // вычисление потенциала и силы, или их обнуление
			if (rij_abs > RCUT) continue;

			vec_normalize(); // подсчет вектора деленного на длину

			// расчет сил взаимодействия
			tx = F * rij_normalize[0]; ty = F * rij_normalize[1]; tz = F * rij_normalize[2];
			Fx[i] += tx; Fx[j] -= tx;
			Fy[i] += ty; Fy[j] -= ty;
			Fz[i] += tz; Fz[j] -= tz;

			Epot += U; // перерасчет потенциальной энергии системы
		}
	}
}

// Расчет ПГУ
void PGU(int index)
{
	if (coordx[index] >= LX) coordx[index] -= LX;
	if (coordx[index] < .0) coordx[index] += LX;
	if (coordy[index] >= LY) coordy[index] -= LY;
	if (coordy[index] < .0) coordy[index] += LY;
	if (coordz[index] >= LZ) coordz[index] -= LZ;
	if (coordz[index] < .0) coordz[index] += LZ;
}

void Verle_scheme() 
{
	Epot = 0; // потенциальная энергия

	// Обновление координат частиц
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		coordx[i] = coordx[i] + vx[i] * STEP + Fx[i] * STEP * STEP / (2 * MASS);
		coordy[i] = coordy[i] + vy[i] * STEP + Fy[i] * STEP * STEP / (2 * MASS);
		coordz[i] = coordz[i] + vz[i] * STEP + Fz[i] * STEP * STEP / (2 * MASS);
		PGU(i);
	}

	// Полушаг обновления скоростей
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		vx[i] += Fx[i] / (2 * MASS) * STEP;
		vy[i] += Fy[i] / (2 * MASS) * STEP;
		vz[i] += Fz[i] / (2 * MASS) * STEP;
	}

	// Перерасчет сил взаимодействия
	double tx, ty, tz;
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		tx = 0; ty = 0; tz = 0;
		// перерасчет сил взаимодействия
		for (int j = 0; j < NUMBERPARTICLES; ++j)
		{
			calc_particle_interaction(i, j); // вычисление потенциала и силы, или их обнуление
			if (rij_abs > RCUT || i == j) continue;

			vec_normalize(); // подсчет вектора деленного на длину

			// расчет сил взаимодействия		
			tx += F * rij_normalize[0];
			ty += F * rij_normalize[1];
			tz += F * rij_normalize[2];

			Epot += U; // перерасчет потенциальной энергии системы
		}
		Fx[i] = tx;
		Fy[i] = ty;
		Fz[i] = tz;
	}

	Epot /= 2; // Избавление от повторяющихся значений

	// конечный шаг расчета скоростей
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		vx[i] += Fx[i] / (2 * MASS) * STEP;
		vy[i] += Fy[i] / (2 * MASS) * STEP;
		vz[i] += Fz[i] / (2 * MASS) * STEP;
	}
}

void Verle_scheme_virtual()
{
	Epot = 0; // потенциальная энергия

	// Обновление координат частиц
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		coordx[i] = coordx[i] + vx[i] * STEP + Fx[i] * STEP * STEP / (2 * MASS);
		coordy[i] = coordy[i] + vy[i] * STEP + Fy[i] * STEP * STEP / (2 * MASS);
		coordz[i] = coordz[i] + vz[i] * STEP + Fz[i] * STEP * STEP / (2 * MASS);
		PGU(i);
	}

	// Полушаг обновления скоростей
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		vx[i] += Fx[i] / (2 * MASS) * STEP;
		vy[i] += Fy[i] / (2 * MASS) * STEP;
		vz[i] += Fz[i] / (2 * MASS) * STEP;
	}

	// Перерасчет сил взаимодействия
	double tx, ty, tz;
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		tx = 0; ty = 0; tz = 0;

		// перерасчет сил взаимодействия
		for (int j = 0; j < NUMBERPARTICLES; ++j)
		{	
			for (int dx = -1; dx < 2; ++dx) {
				for (int dy = -1; dy < 2; ++dy) {
					for (int dz = -1; dz < 2; ++dz) {
						if (i == j && dx == 0 && dy == 0 && dz == 0) continue;
						double* vec_tmp = new double[3] {coordx[i] - (coordx[j] + dx * LX), coordy[i] - (coordy[j] + dy * LY), coordz[i] - (coordz[j] + dz * LZ)};
						double vec_tmp_abs = sqrt(vec_euklid(vec_tmp));
						std::pair<double, double>* u_and_f = calc_particle_interaction(vec_tmp); // расчет потенциала им силы взаимодействия
						double tmp_F = u_and_f->second;
						if (vec_tmp_abs > RCUT) {
							delete u_and_f;
							delete[] vec_tmp;
							continue;
						}
						tx += tmp_F * vec_tmp[0] / vec_tmp_abs;
						ty += tmp_F * vec_tmp[1] / vec_tmp_abs;
						tz += tmp_F * vec_tmp[2] / vec_tmp_abs;

						Epot += u_and_f->first; // перерасчет потенциальной энергии системы

						delete u_and_f;
						delete[] vec_tmp;
					}

				}

			}
			
		}
		Fx[i] = tx;
		Fy[i] = ty;
		Fz[i] = tz;
	}

	Epot /= 2; // Избавление от повторяющихся значений

	// конечный шаг расчета скоростей
	for (int i = 0; i < NUMBERPARTICLES; ++i)
	{
		vx[i] += Fx[i] / (2 * MASS) * STEP;
		vy[i] += Fy[i] / (2 * MASS) * STEP;
		vz[i] += Fz[i] / (2 * MASS) * STEP;
	}
}

// Подсчет средней скорости (скорость центра масс)
double* calc_vm()
{
	double sumX = .0, sumY = .0, sumZ = .0;
	for (int i = 0; i < NUMBERPARTICLES; i++)
	{
		sumX += vx[i];
		sumY += vy[i];
		sumZ += vz[i];
	}	
	return new double[3] {sumX /= NUMBERPARTICLES, sumY /= NUMBERPARTICLES, sumZ /= NUMBERPARTICLES};
}

void calc_Ekin_Epot_P()
{
	P = 0; // давление нуль
	double* vi = new double[3] { 0., 0., 0. };
	double* tmp = new double[3] { 0., 0., 0. }; //  r_ij * F_ij (для давления)

	double sum_vi_kin = 0; // сумма квадратов координат скоростей для кинетической энергии
	double sum_vi_term = 0; // сумма квадратов координат скоростей для тепловой энергии
	double sum_for_P = 0; // сумма r_ij * F_ij (для давления) 
	double* vm = calc_vm(); // скорость центра масс 

	for (int i = 0; i < NUMBERPARTICLES; ++i) {
		vi[0] = vx[i]; vi[1] = vy[i]; vi[2] = vz[i];
		sum_vi_kin += vec_euklid(vi);
		vi[0] -= vm[0]; vi[1] -= vm[1]; vi[2] -= vm[2];
		sum_vi_term += vec_euklid(vi);

		// подсчет правой части уравнения для давления
		for (int j = i + 1; j < NUMBERPARTICLES; ++j) {
			vec_rij(i, j);
			vec_absolute();
			vec_normalize();
			if (rij_abs <= RCUT) {
				calc_particle_interaction(i, j);
				tmp[0] = F * rij_normalize[0] * rij[0];
				tmp[1] = F * rij_normalize[1] * rij[1];
				tmp[2] = F * rij_normalize[2] * rij[2];
				sum_for_P += tmp[0] + tmp[1] + tmp[2];
			}
		}
	}
	Ekin = MASS * sum_vi_kin / 2;
	Eterm = MASS * sum_vi_term / 2;
	P = (MASS * sum_vi_term + sum_for_P) / (3 * VOLUME);

	delete[] vi;
	delete[] tmp;
	delete[] vm;
}

void upd_E() { E = Ekin + Epot; }
void upd_Eint() { Eint = Eterm + Epot; }
void upd_T() { T = 2 * Eterm / (3 * NUMBERPARTICLES * K_B); }

// 10 лабораторная работа
void task_10()
{
	using namespace std;
	initial_condition_two_particles();
	initial_vectors_F();

	// работа с файлами
	ofstream stream;

	stream.open("C:\\Users\\Dima\\source\\repos\\Molecular Dinamic\\Gimazetdinov_MD_10.txt");
	stream << fixed;
	stream.precision(8); // установка глобальной точности



	for (int i = 0; i < NSTEPS; ++i)
	{

		if (stream.is_open() && (i < 1000))
		{
			vec_rij(0, 1);
			vec_normalize();
			vec_absolute();

			stream << "Step = " << i << endl;
			stream << "r1 = (rx1;ry1;rz1) = (" << coordx[0] << ";" << coordy[0] << ";" << coordz[0] << ")" << endl;
			stream << "r2 = (rx2;ry2;rz2) = (" << coordx[1] << ";" << coordy[1] << ";" << coordz[1] << ")" << endl;
			stream << "r12_abs = " << rij_abs << endl;

			calc_particle_interaction(0, 1);
			stream << "U12 = " << U << endl;
			stream << "F12 = " << F << endl;
			stream << "F1 = (Fx1;Fy1;Fz1) = (" << Fx[0] << ";" << Fy[0] << ";" << Fz[0] << ")" << endl;
			stream << "v1=(vx1;vy1;vz1) = (" << vx[0] << ";" << vy[0] << ";" << vz[0] << ")" << std::endl;
			stream << "v2=(vx2;vy2;vz2) = (" << vx[1] << ";" << vy[1] << ";" << vz[1] << ")" << std::endl << std::endl;
		}

		// перерасчет параметров c ПГУ
		Verle_scheme();
	}

	// закрытие потока записи и высвобождение памяти
	stream.close();
}

void task_11()
{
	using namespace std;

	initial_condition_two_particles();
	initial_vectors_F();
	calc_Ekin_Epot_P(); // расчет кинетической и тепловой
	upd_E(); // расчет полной энергии
	upd_Eint(); // расчет полной внутренней энергии
	upd_T(); // расчет температуры системы

	ofstream stream;
	ofstream plot;
	stream.open("C:\\Users\\Dima\\source\\repos\\Molecular Dinamic\\Gimazetdinov_MD_11.txt");
	stream << fixed;
	stream.precision(8);

	plot.open("C:\\Users\\Dima\\source\\repos\\MD_Full_Rework\\MD_Full_Rework\\Dots.txt");
	plot << fixed;
	plot.precision(8);

	for (int i = 0; i < NSTEPS; i++)
	{
		if (plot.is_open()) plot << Ekin << ',' << Epot << ',' << E << ',' << T << ',' << P << ',';
		if (stream.is_open())
		{
			stream << "Step=" << i << endl;
			stream << "Ekin=" << Ekin << endl;
			stream << "Eterm=" << Eterm << endl;
			stream << "Epot=" << Epot << endl;
			stream << "Eint=" << Eint << endl;
			stream << "E=" << E << endl;
			stream << "T=" << T << endl;
			stream << "P=" << P << endl << endl;
		}
		// Перерасчет параметров
		Verle_scheme();
		
		
		calc_Ekin_Epot_P(); // расчет кинетической и тепловой
		upd_E(); // расчет полной энергии
		upd_Eint(); // расчет полной внутренней энергии
		upd_T(); // расчет температуры системы
	}

	stream.close();
	plot.close();
}

void task_12()
{
	using namespace std;

	initial_condition_two_particles();
	initial_vectors_F_virt();
	calc_Ekin_Epot_P(); // расчет кинетической и тепловой
	upd_E(); // расчет полной энергии
	upd_Eint(); // расчет полной внутренней энергии
	upd_T(); // расчет температуры системы

	ofstream stream;
	ofstream plot1, plot2;
	stream.open("C:\\Users\\Dima\\source\\repos\\Molecular Dinamic\\Gimazetdinov_MD_12.txt");
	stream << fixed;
	stream.precision(8);

	plot1.open("C:\\Users\\Dima\\source\\repos\\Molecular Dinamic\\p1.txt");
	plot1 << fixed;
	plot1.precision(8);

	plot2.open("C:\\Users\\Dima\\source\\repos\\Molecular Dinamic\\p2.txt");
	plot2 << fixed;
	plot2.precision(8);

	for (int i = 0; i < NSTEPS; i++)
	{
		if (plot1.is_open()) plot1 << coordx[0] << endl << coordy[0] << endl;
		if (plot2.is_open()) plot2 << coordx[1] << endl << coordy[1] << endl;
		if (stream.is_open())
		{
			vec_rij(0, 1);
			vec_normalize();
			vec_absolute();

			stream << "Step = " << i << endl;
			stream << "r1 = (rx1;ry1;rz1) = (" << coordx[0] << ";" << coordy[0] << ";" << coordz[0] << ")" << endl;
			stream << "r2 = (rx2;ry2;rz2) = (" << coordx[1] << ";" << coordy[1] << ";" << coordz[1] << ")" << endl;
			stream << "r12_abs = " << rij_abs << endl;

			calc_particle_interaction(0, 1);
			stream << "U12 = " << U << endl;
			stream << "F12 = " << F << endl;
			stream << "F1 = (Fx1;Fy1;Fz1) = (" << Fx[0] << ";" << Fy[0] << ";" << Fz[0] << ")" << endl;
			stream << "v1=(vx1;vy1;vz1) = (" << vx[0] << ";" << vy[0] << ";" << vz[0] << ")" << std::endl;
			stream << "v2=(vx2;vy2;vz2) = (" << vx[1] << ";" << vy[1] << ";" << vz[1] << ")" << std::endl << std::endl;

			stream << "Step=" << i << endl;
			stream << "Ekin=" << Ekin << endl;
			stream << "Eterm=" << Eterm << endl;
			stream << "Epot=" << Epot << endl;
			stream << "Eint=" << Eint << endl;
			stream << "E=" << E << endl;
			stream << "T=" << T << endl;
			stream << "P=" << P << endl << endl;
		}
		// Перерасчет параметров
		Verle_scheme_virtual();

		// где-то ниже утечка памяти
		calc_Ekin_Epot_P(); // расчет кинетической и тепловой
		upd_E(); // расчет полной энергии
		upd_Eint(); // расчет полной внутренней энергии
		upd_T(); // расчет температуры системы
	}

	stream.close();
	plot1.close();
	plot2.close();
}

int main() {
	memory_allocation();
	task_12();
	memory_free();

	delete[] rij;
	delete[] rij_normalize;

	return 0;
}