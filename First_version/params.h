#pragma once

const int NSTEPS = 1000; // ����� �����
const int LASTSTEP = NSTEPS - 1; // ��������� ���

// ����� ������������ ����� (����������) �� ���� ��������� 
const int NUMCRIST_X = 2; // ����� ������������ ����� �� ��� X
const int NUMCRIST_Y = 2; // ����� ������������ ����� �� ��� Y
const int NUMCRIST_Z = 2; // ����� ������������ ����� �� ��� Z


// ����� ������ ��� ����������� ������� � ������ ��� (������������� ��������� �������)
// const int NUMBERPARTICLES = NUMCRIST_X * NUMCRIST_Y * NUMCRIST_Z;
const int NUMBERPARTICLES = 2; //����� ������ 
const double STEP = 0.002; // ��� �������������� ���������� �����