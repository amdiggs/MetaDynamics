#ifndef AMDmath_hpp
#define AMDmath_hpp


#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <thread>
#include <chrono>
#include <OpenCL/opencl.h>
#define twoPI 6.283185307
#define PI 3.141592654
#define halfPI 1.5707963

#define U_2_G 1.66e+24


namespace AMD {
struct Vec3;

int Min(int a, int b);
float Min(float a, float b);
int Max(int a, int b);
float Max(float a, float b);

float AMU_to_G(float u);

float LERP(float val, float max, float min);
float Round(float val, int num_decimal);

int sgn(int a);
float sgn(float a);


struct Vec3{
    float x,y,z;
    Vec3();
    Vec3(float s);
    Vec3(float e_x, float e_y, float e_z);
    
    Vec3 add(const Vec3& other) const;
    float& operator[](const int index);
    float* get();
    
    float dot(const Vec3& other) const;
    Vec3 cross(const Vec3& other) const;
    float len() const;
    
    
    Vec3 operator+(const Vec3& other) const;
    Vec3 operator-(const Vec3& other) const;
    Vec3& operator=(const Vec3& other);
    operator cl_float3() const;
    Vec3 operator*(const Vec3& other) const;
    Vec3 operator/(float div) const;
    Vec3 operator+=(const Vec3& other);
    Vec3 operator*=(float scale);
    
    bool operator==(const Vec3& other) const;
    bool Is_Parallel(const Vec3& other) const;
    void Reset();
    void print();
    void Normalize();
    void Vround(int decimals);
    
    
};

Vec3 operator*(const float scale, const Vec3& other);
Vec3 operator*(const Vec3& other, const float scale);
Vec3 Normalize(const Vec3& vec);
Vec3 Round(const Vec3& vec, int decimals);
float Get_angle(const Vec3& A, const Vec3& B);

float Distance(const Vec3& A, const Vec3& B);

// ############# Matrix ########################

struct Mat2
{
    float m[2][2];
    Mat2();
    Mat2(float a,float b,float c,float d);
    float Det();
};

struct Mat3{
    const unsigned int n =3;
    float m[3][3];
    Mat3();
    
    
    float* operator[](const int& index);
    Mat3& operator=(const Mat3& other);
    Mat3 add(const Mat3& other) const;
    Mat3 operator+(const Mat3& other) const;
    
    Mat3 multiply(const Mat3& other) const;
    Mat3 operator*(const Mat3& other) const;
    Vec3 multiply(Vec3& other) const;
    Vec3 operator*(Vec3& other) const;
    
    void assign_col(int col_idx, Vec3 col);
    void assign_row(int row_idx, Vec3 row);
    
    void Rotate(Vec3 ang);
    void Scale(float);
    void Scale(Vec3);
    void Transpose();
    float Det() const;
    Mat2 Adj(int i, int j) const;
    void print();
};

Mat3 Inverse(const Mat3& other);





} //end of namespace
#endif
