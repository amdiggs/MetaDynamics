
#include "MyVec.hpp"
#include <cstdio>
using namespace AMD;



int AMD::Min(int a, int b){
    return (a<b) ? a : b;
}

float AMD::Min(float a, float b){
    return (a<b) ? a : b;
}

int AMD::Max(int a, int b){
    return (a>b) ? a : b;
}

float AMD::Max(float a, float b){
    return (a>b) ? a : b;
}


float AMD::AMU_to_G(float u){
    return u/U_2_G;
}
float AMD::Round(float val, int num_decimal){
    float mod = 10.0;
    for( int i = 0; i<num_decimal; i++){
        mod*=10.0;
    }
    int temp = (int)(val*mod);
    return (float)temp / mod;
}

int AMD::sgn(int a)
{
    return (a > 0)? 1 : -1;
}

float AMD::sgn(float a)
{
    return (a>0.0)? 1.0 : -1.0;
}
//           __     ___  ___
//  \    /  |      /        \
//   \  /   |--   |       __/
//    \/    |__    \___     \
//                        __/
Vec3::Vec3()
:x(0.0), y(0.0), z(0.0)
{}


Vec3::Vec3(float s)
:x(s) , y(s) , z(s)
{}

Vec3::Vec3(float e_x, float e_y, float e_z)
:x(e_x), y(e_y),z(e_z)
{}


float* Vec3::get(){
    return &x;
}


float& Vec3::operator[](const int index){
    return get()[index];
}

float Vec3::dot(const AMD::Vec3 &other) const {
    return x*other.x + y*other.y + z*other.z;
}

AMD::Vec3 Vec3::cross(const AMD::Vec3 &other) const {
    float _x = y*other.z - z*other.y;
    float _y = z*other.x - x*other.z;
    float _z = x*other.y - y*other.x;
    return Vec3(_x, _y, _z);
}

float Vec3::len() const{
    return sqrt(this->dot(*(this)));
}

AMD::Vec3 &Vec3::operator=(const AMD::Vec3 &other) {
    if (&other == this){
        return  *this;
    }
    else {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
    }
    return *this;
}




AMD::Vec3::operator cl_float3() const {
    cl_float3 ret;
    ret.s[0] = x;
    ret.s[1] = y;
    ret.s[2] = z;
    return ret;
}


AMD::Vec3 Vec3::operator*(const AMD::Vec3 &other) const{
    return Vec3(x*other.x,y*other.y, z*other.z);
}



AMD::Vec3 Vec3::operator/(float div) const{
    return Vec3(x/div,y/div,z/div);
}


AMD::Vec3 Vec3::operator-(const Vec3& other) const{
    float _x = this->x - other.x;
    float _y = this->y - other.y;
    float _z = this->z - other.z;
    return Vec3(_x,_y,_z);
    
}




AMD::Vec3 Vec3::operator+(const Vec3& other) const{
    float _x = this->x + other.x;
    float _y = this->y + other.y;
    float _z = this->z + other.z;
    return Vec3(_x,_y,_z);
    
}




AMD::Vec3 Vec3::operator+=(const Vec3& other){
    this->x = this->x + other.x;
    this->y = this->y + other.y;
    this->z = this->z + other.z;
    
    return *this;
}

AMD::Vec3 Vec3::operator*=(float scale){
    *this = *this * scale;
    return *this;
    
}

bool Vec3::operator==(const Vec3& other) const{
    return (this->x == other.x && this->y == other.y && this->z == other.z);
}


bool Vec3::Is_Parallel(const Vec3& other) const{
    Vec3 rev = other*(-1.0);
    return (*this == other || *this == rev);
}

void AMD::Vec3::print(){
    std::cout << x << " " << y << " " << z << std::endl;
}


void AMD::Vec3::Reset(){
    x = 0.0;
    y = 0.0;
    z = 0.0;
}

void AMD::Vec3::Normalize(){
    float N = this->len();
    *this = (*this)/N;
}




void AMD::Vec3::Vround(int decimals){
    x = Round(x, decimals);
    y = Round(y, decimals);
    z = Round(z, decimals);
}

AMD::Vec3 AMD::operator*(const float scale, const Vec3& other){
    return Vec3(other.x*scale,other.y*scale, other.z*scale);
}

AMD::Vec3 AMD::operator*(const Vec3& other, const float scale){
    return Vec3(other.x*scale,other.y*scale, other.z*scale);
}

AMD::Vec3 AMD::Normalize(const AMD::Vec3& vec){
    float x = vec.x;
    float y = vec.y;
    float z = vec.z;
    float N = sqrt(x*x + y*y + z*z);
    return AMD::Vec3(x/N, y/N, z/N);
}




float AMD::Get_angle(const Vec3& A, const Vec3& B){
    if(A.len() <= 0.00001 || B.len() <= 0.00001){
        return 1.57079;
    }
    else{
        float c_th = A.dot(B)/(A.len()*B.len());
        return acos(c_th);
    }
}



float AMD::Distance(const Vec3& A, const Vec3& B){
    AMD::Vec3 temp = A - B;
    return temp.len();
}

AMD::Vec3 AMD::Round(const Vec3& vec, int decimals){
    return AMD::Vec3(Round(vec.x, decimals),Round(vec.y, decimals),Round(vec.z, decimals));
}


//#################################################################################
//#################################################################################
//#################################################################################

//==============MATRICIES==========================================


Mat2::Mat2()
{
    m[0][0]=0.;
    m[0][1]=0.;
    m[1][0]=0.;
    m[1][1]=0.;
}


Mat2::Mat2(float a, float b, float c, float d)
{
    m[0][0]=a;
    m[0][1]=b;
    m[1][0]=c;
    m[1][1]=d;
}

float Mat2::Det()
{
    return m[0][0]*m[1][1] - m[1][0]*m[0][1];
}

Mat3::Mat3() {
    for (int i = 0; i< 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j){
                m[i][j] = 1.0;
            }
            else{m[i][j] = 0.0;}
        }
    }
}

float *Mat3::operator[](const int &index) {
    return m[index];
}

AMD::Mat3 &Mat3::operator=(const AMD::Mat3 &other) {
    if(this == &other){ return  *this;}
    else{
        for (int i = 0; i< 3; i++) {
            for (int j = 0; j < 3; j++) {
            this->m[i][j] = other.m[i][j];
            }

        }
    
        return *this;
    }
}

AMD::Mat3 Mat3::add(const AMD::Mat3 &other) const {
    Mat3 temp;
    for (int i =0; i<n; i++){
        for (int j = 0; j < n; j++) {
            temp.m[i][j] = m[i][j] + other.m[i][j];
        }
    }
    return temp;
}

AMD::Mat3 Mat3::operator+(const AMD::Mat3 &other) const {
    return this->add(other);
}

AMD::Mat3 Mat3::multiply(const AMD::Mat3 &other) const {
    Mat3 temp;
    float el;
    for (int col = 0; col<n; col++){
        for (int i =0; i<n; i++){
            el =0.0;
            for (int j = 0; j < n; j++) {
                el += other.m[i][j] * m[j][col];
            }
            temp.m[i][col] = el;
        }
    }
    return temp;
}

AMD::Mat3 Mat3::operator*(const AMD::Mat3 &other) const {
    return this->multiply(other);
}

AMD::Vec3 Mat3::multiply(AMD::Vec3 &other) const {
    Vec3 temp;
    for (int row = 0; row<n; row++){
        for (int col = 0; col<n;col++){
            temp[row] += m[col][row]*other[col];
        }
    }
    return temp;
}

AMD::Vec3 Mat3::operator*(AMD::Vec3 &other) const {
    return this->multiply(other);
}

void Mat3::assign_col(int col_idx, AMD::Vec3 col) {
    for (int i = 0; i<n; i++){
        m[i][col_idx] = col[i];
    }
}

void Mat3::assign_row(int row_idx, AMD::Vec3 row) {
    for (int i = 0; i<n; i++){
        m[row_idx][i] = row[i];
    }
}

void Mat3::Rotate(AMD::Vec3 ang) {
    float a = ang.x; float b = ang.y; float c = ang.z;
    AMD::Mat3 rot;
    AMD::Vec3 c0 (cos(a)*cos(b),  cos(a)*sin(b)*sin(c) - sin(a)*cos(c),  cos(a)*sin(b)*cos(c)+sin(a)*sin(c));
    AMD::Vec3 c1(sin(a)*cos(b), sin(a)*sin(b)*sin(c) + cos(a)*cos(c), sin(a)*sin(b)*cos(c)-cos(a)*sin(c));
    AMD::Vec3 c2(-sin(b), cos(b)*sin(c), cos(b)*cos(c));
    
    rot.assign_col(0, c0);
    rot.assign_col(1, c1);
    rot.assign_col(2, c2);
    *this = rot * (*this);
    
    return;
}

void Mat3::Scale(float scale) {
    Mat3 temp;
    temp[0][0] = scale;
    temp[1][1] = scale;
    temp[2][2] = scale;
    
    *this = temp * (*this);
    
    return;
}

void Mat3::Scale(AMD::Vec3 vec) {
    
    m[0][0] = m[0][0] * vec.x;
    m[1][1] = m[1][1] * vec.y;
    m[2][2] = m[2][2] * vec.z;
    
    
    return;
}

void Mat3::Transpose(){
    Mat3 temp = *this;
    
    for (int i = 0; i< n; i++) {
        for (int j = 0; j < n; j++) {
                m[i][j] = temp.m[j][i];
        }
    }
}


AMD::Mat2 Mat3::Adj(int r, int c) const
{
    AMD::Vec3 rows[2];
    int count = 0;
    for(int i = 0; i<3;i++){
	if(i == r){continue;}
	else{
	    rows[count] = AMD::Vec3(m[i][0],m[i][1],m[i][2]);
	    count++;
	}
    }
    float r1[2], r2[2];
    count = 0;
    for(int i = 0; i<3;i++){
	if(i == c){continue;}
	else{
	    r1[count] = rows[0][i];
	    r2[count] = rows[1][i];
	    count++;
	}
    }
    return Mat2(r1[0],r1[1],r2[0],r2[1]);
}

float Mat3::Det() const
{
    float ret = 0.;
    float sgn = 1.0;
    for(int i = 0; i<3; i++){
	Mat2 adj = Adj(0,i);
	ret += sgn*m[0][i]*adj.Det();
	sgn*=-1.0;
    }
    return ret;
}


void Mat3::print() {
    std::cout <<  std::endl;
    for (int i =0; i<n; i++){
        for (int j = 0; j < n; j++) {
        std::cout << m[i][j] << ", ";
        }
       std::cout <<  std::endl;
    
    }
    std::cout <<  std::endl;
}



AMD::Mat3 AMD::Inverse(const Mat3& other)
{
    AMD::Mat3 ret;
    float sgn = 1.0;
    float det = other.Det();
    if(abs(det) <= 0.000001){printf("Matrix is not invertable!!!/n"); exit(10);}
    for(int i = 0; i<3; i++){
	    for(int j = 0; j<3; j++){
	    AMD::Mat2 adj = other.Adj(i,j);
	    float el = sgn*adj.Det();
	    ret[j][i] = el/det;
	    }
    }

    return ret;
}


