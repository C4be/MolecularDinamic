#pragma once
#include <cmath>
#include "params.h"

const double MASS = 66.335;  // ����� ���� ������� (���� ������)
const double K_B = 1.380648528;  // ���������� ���������

// ��������� ���������� �.-��. ��� ������
const double EPS = 1.712;
const double SIGMA = 0.3418; // �������� ����� ��������������
const double RCUT = 2.5 * SIGMA; // ������ ��������� ����������
const double RCUT2 = RCUT * RCUT;
const double UCUT = 4 * EPS * (pow((SIGMA / RCUT), 12) - pow((SIGMA / RCUT), 6)); // ��������� ���������
//const double UCUT = 0;
const double ACRIST = 1.0; // ����� ����� ��. ������, ������� �� ������������������ ��������� � ������

// ������� ������� �� ���� ���������
const double LX = NUMCRIST_X * ACRIST, LY = NUMCRIST_Y * ACRIST, LZ = NUMCRIST_Z * ACRIST;
const double VOLUME = LX * LY * LZ; // ����� �������
