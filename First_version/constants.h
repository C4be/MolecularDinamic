#pragma once
#include <cmath>
#include "params.h"

const double MASS = 66.335;  // масса одно частицы (атом Аргона)
const double K_B = 1.380648528;  // Постоянная Больцмана

// параметры потенциала Л.-Дж. для аргона
const double EPS = 1.712;
const double SIGMA = 0.3418; // параметр длины взаимодействия
const double RCUT = 2.5 * SIGMA; // радиус обрезания потенциала
const double RCUT2 = RCUT * RCUT;
const double UCUT = 4 * EPS * (pow((SIGMA / RCUT), 12) - pow((SIGMA / RCUT), 6)); // потенциал обрезания
//const double UCUT = 0;
const double ACRIST = 1.0; // длина ребра эл. ячейки, зависит от термодинамического состояния и модели

// размеры системы по осям координат
const double LX = NUMCRIST_X * ACRIST, LY = NUMCRIST_Y * ACRIST, LZ = NUMCRIST_Z * ACRIST;
const double VOLUME = LX * LY * LZ; // объем системы
