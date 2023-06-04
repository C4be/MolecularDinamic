#pragma once

const int NSTEPS = 1000; // число шагов
const int LASTSTEP = NSTEPS - 1; // последний шаг

// число элементарных €чеек (кристаллов) по ос€м координат 
const int NUMCRIST_X = 2; // число элементарных €чеек по оси X
const int NUMCRIST_Y = 2; // число элементарных €чеек по оси Y
const int NUMCRIST_Z = 2; // число элементарных €чеек по оси Z


// число частиц дл€ примитивной решетки с учетом ѕ√” (периодические граничные услови€)
// const int NUMBERPARTICLES = NUMCRIST_X * NUMCRIST_Y * NUMCRIST_Z;
const int NUMBERPARTICLES = 2; //число частиц 
const double STEP = 0.002; // шаг интегрировани€ разностной схемы