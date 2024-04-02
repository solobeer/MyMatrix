#ifndef MYMATRIX_H
#define MYMATRIX_H
#include <iostream>  
#include <cstdlib>  
#include <cmath>
#include <vector>
#include <typeinfo>
#include <omp.h>
//#include <intrin.h>
//#include <immintrin.h>

namespace MyMat {
	template  <class T>
	class MyMatrix
	{
	private:
		int cols;
		int rows;
		int size;
		T* data;
		bool is_singular = false;
	public:
		MyMatrix();//Ĭ�Ϲ��캯��
		MyMatrix(int rows_, int cols_);
		MyMatrix(int rows_, int cols_, T val);//���������������,Ԫ�����ֵval
		MyMatrix(T* buffer, int rows_, int cols_, bool isRowFirst = true);
		MyMatrix(const MyMatrix& MyMatrix);                     //���Ѵ��ڵľ�����������  
		inline int getCols() const { return cols; };                //��ȡ����  
		inline int getRows() const { return rows; };                //��ȡ����  
		inline int getSize() const { return rows * cols; };         //��ȡ�����С
		bool isSingular() { inverse(); return is_singular; };//�жϾ����Ƿ�����
		void setValue(int rows_, int cols_, T value);
		T getValue(int rows_, int cols_);

		T& operator()(int row, int col);					//���Ų��������أ����ڻ�ȡ�����row�е�col��Ԫ�� 

		bool LUP_Descomposition(MyMatrix<double>& A, MyMatrix<double>& L, MyMatrix<double>& U, MyMatrix<double>& P);//LUP�ֽ�
		MyMatrix<double> LUP_Solve(MyMatrix<double>& L, MyMatrix<double>& U, MyMatrix<double>& P, MyMatrix<double>& b);//LUP������Է���

		MyMatrix<double>  inverse();                                  //LU�ֽ�������� 
		MyMatrix<T>& operator=(const MyMatrix& MyMatrix);            //����ĸ�ֵ����������  
		MyMatrix<T> MyMatrix<T>::mul(const MyMatrix<T>& MyMatrix1, const MyMatrix<T>& MyMatrix2);			//�������

		~MyMatrix() { delete[]data; };
	public:
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator+(const MyMatrix<ElemType>& Matrix1, const MyMatrix<ElemType>& Matrix2);//���������ӷ�+����������
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator-(const MyMatrix<ElemType>& Matrix1, const MyMatrix<ElemType>& Matrix2);//������������-����������
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator*(const MyMatrix<ElemType>& Matrix1, const MyMatrix<ElemType>& Matrix2);//���������˷�*����������

	
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator+(const MyMatrix<ElemType>& Matrix, const ElemType val);      //���������ӷ�+���������أ�1��
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator+(const ElemType val, const MyMatrix<ElemType>& Matrix);      //���������ӷ�+���������أ�2��
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator-(const MyMatrix<ElemType>& Matrix, const ElemType val);      //������������-���������أ�1��
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator-(const ElemType val, const MyMatrix<ElemType>& Matrix);      //������������-���������أ�2��
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator*(const MyMatrix<ElemType>& Matrix, const ElemType val);      //���������˷�*���������أ�1��
		template  <class ElemType>
		friend MyMatrix<ElemType>  operator*(const ElemType val, const MyMatrix<ElemType>& Matrix);      //���������˷�*���������أ�2��

		template  <class ElemType>
		friend std::ostream& operator<<(std::ostream &os, const MyMatrix<ElemType>& MyMatrix); //��������������
	};

	//Ĭ�Ϲ��캯��
	template  <class T>
	MyMatrix<T>::MyMatrix()
	{
		cols = 0;
		rows = 0;
		size = 0;
		data = nullptr;
	}
	template <class T>
	MyMatrix<T>::MyMatrix(int rows_, int cols_) {
		cols = cols_;
		rows = rows_;
		size = cols*rows;
		data = new T[size];
		for (int i = 0; i < size; i++)
			data[i] = T(0);
	}

	//���캯���������й���
	template  <class T>
	MyMatrix<T>::MyMatrix(int rows_, int cols_, T val)
	{
		cols = cols_;
		rows = rows_;
		size = cols*rows;
		data = new T[size];
		if (data == nullptr) cerr << "allocate error" << endl;
		for (int i = 0; i < size; i++)
			data[i] = val;
	}

	template <class T>
	MyMatrix<T>::MyMatrix(T* buffer, int rows_, int cols_, bool isRowFirst = true) {
		cols = cols_;
		rows = rows_;
		size = cols * rows;
		data = new T[size];
		if (isRowFirst == false) {
			/*for (int i = 0; i < size; i++) {
			data[i] = Array[i];
			}*/
		}
		else {
			for (int i = 0; i < size; i++) {
				data[i] = buffer[i];
			}
		}
	}

	//���캯�������๹��  
	template  <class T>
	MyMatrix<T>::MyMatrix(const MyMatrix& MyMatrix)
	{
		cols = MyMatrix.cols;
		rows = MyMatrix.rows;
		size = cols*rows;
		data = new T[size];
		for (int i = 0; i < size; i++)
			data[i] = MyMatrix.data[i];
	}
	template <class T>
	void MyMatrix<T>::setValue(int rows_, int cols_, T value) {
		if (0 <= rows_ && rows_ < rows && 0 <= cols_ && cols_ < cols)
			data[rows_ * cols + cols_] = value;
		else
			std::cerr << "access out of range" << endl;
	}

	template <class T>
	T MyMatrix<T>::getValue(int rows_, int cols_) {
		if (0 <= rows_ && rows_ < rows && 0 <= cols_ && cols_ < cols)
			return data[rows_ * cols + cols_];
		else
		{
			std::cerr << "access out of range" << endl;
			return 0;
		}
	}


	//----------------------------------------------------------------------------------- 
	template  <class T>
	T& MyMatrix<T>::operator()(int row, int col)
	{
		return data[col * rows + row];
	}

	//�������  
	template  <class ElemType>
	std::ostream& operator<<(std::ostream& os, const MyMatrix<ElemType>& MyMatrix)
	{
		for (int i = 0; i < MyMatrix.rows; i++)
		{
			for (int j = 0; j < MyMatrix.cols; j++)
			{
				os << MyMatrix.data[i * MyMatrix.cols + j];
				if (j != MyMatrix.cols - 1)
					os << ",";
			}
			os << ";" << std::endl;
		}
		return os;
	}

	//������ = ����
	template  <class T>
	MyMatrix<T>& MyMatrix<T>::operator=(const MyMatrix& MyMatrix)
	{
		//һ��Ҫ������ʱ���󣬷�ֹͬһ����ֵ�����·����ڴ��ͬʱ���Ѿ����˾���ı䣡������������  
		if (this == &MyMatrix)
		{
			return *this;
		}
		cols = MyMatrix.cols;
		rows = MyMatrix.rows;
		size = cols*rows;
		delete[] data;
		data = new T[size];
		for (int i = 0; i < size; i++)
			data[i] = MyMatrix.data[i];
		return *this;
	}

	//������ + ����(�����)  
	template  <class ElemType>
	MyMatrix<ElemType>  operator+(const MyMatrix<ElemType>& Matrix1, const MyMatrix<ElemType>& Matrix2)
	{
		MyMatrix<ElemType> res(Matrix1.rows, Matrix1.cols);
		if (Matrix1.cols != Matrix2.cols || Matrix1.rows != Matrix2.rows)
		{
			std::cerr << "Error:The number of rows or columns is not equal(MyMatrix::operator+)" << std::endl;
		}
		else if (typeid(Matrix1.data).name() != typeid(Matrix2.data).name())
		{
			std::cerr << "Error:Different types of Matrix data(MyMatrix::operator+)" << std::endl;
		}
		else
		{
			for (int i = 0; i < Matrix1.size; i++)
				res.data[i] = Matrix1.data[i] + Matrix2.data[i];
		}
		return res;
	}

	//������ - ����(�����)  
	template  <class ElemType>
	MyMatrix<ElemType>  operator-(const MyMatrix<ElemType>& MyMatrix1, const MyMatrix<ElemType>& MyMatrix2)
	{
		MyMatrix<ElemType> res(MyMatrix1.rows,MyMatrix1.cols);
		if (MyMatrix1.cols != MyMatrix2.cols || MyMatrix1.rows != MyMatrix2.rows)
		{
			std::cerr << "Error:The number of rows or columns is not equal(MyMatrix::operator-)" << std::endl;
		}
		else if (typeid(MyMatrix1.data).name() != typeid(MyMatrix2.data).name())
		{
			std::cerr << "Error:Different types of Matrix data(MyMatrix::operator-)" << std::endl;
		}
		else
		{
			for (int i = 0; i < MyMatrix1.size; i++)
				res.data[i] = MyMatrix1.data[i] - MyMatrix2.data[i];
		}
		return res;
	}
	//������ * ����(�����)  

	template  <class T>
	MyMatrix<T> MyMatrix<T>::mul(const MyMatrix<T>& MyMatrix1, const MyMatrix<T>& MyMatrix2)
	{
		MyMatrix<T> res(MyMatrix1.rows, MyMatrix2.cols, 0.0);
		if (MyMatrix1.cols != MyMatrix2.rows) {
			std::cerr << "Error:The cols of MyMatrix1 is not equal to the rows of MyMatrix2, and can not be multiplied.(MyMatrix::operator*)" << std::endl;
			return res;
		}
		else if (typeid(MyMatrix1.data).name() != typeid(MyMatrix2.data).name()) {
			std::cerr << "Error:Different types of Matrix data(MyMatrix::operator*)" << std::endl;
			return res;
		}
		else
		{
			if (MyMatrix1.rows >= 65 && MyMatrix1.cols >= 65 && MyMatrix2.cols >= 65) {
				int num = 0;
				if (MyMatrix1.rows < 2 * omp_get_num_procs()) num = MyMatrix1.rows;
				else
					num = 64;
				#pragma omp parallel for num_threads(num)
				for (int i = 0; i < MyMatrix1.rows; i++) {
					//printf("i = %d, I am Thread %d\n", i, omp_get_thread_num());
					for (int k = 0; k < MyMatrix1.cols; k++)
					{
						for (int j = 0; j < MyMatrix2.cols; j++)
						{
							res.data[i * res.cols + j] += MyMatrix1.data[MyMatrix1.cols * i + k] * MyMatrix2.data[k * MyMatrix2.cols + j];
						}
					}
				}
			}
			else {
				for (int i = 0; i < MyMatrix1.rows; i++) {
					/*printf("i = %d, I am Thread %d\n", i, omp_get_thread_num());*/
					for (int k = 0; k < MyMatrix1.cols; k++)
					{
						for (int j = 0; j < MyMatrix2.cols; j++)
						{
							res.data[i * res.cols + j] += MyMatrix1.data[MyMatrix1.cols * i + k] * MyMatrix2.data[k * MyMatrix2.cols + j];
						}
					}
				}
			}
			return res;
		}
	}

	template  <class ElemType>
	MyMatrix<ElemType>  operator*(const MyMatrix<ElemType>& MyMatrix1, const MyMatrix<ElemType>& MyMatrix2)
	{
		MyMatrix<ElemType> res(MyMatrix1.rows, MyMatrix2.cols);
		if (MyMatrix1.cols != MyMatrix2.rows) {
			std::cerr << "Error:The cols of MyMatrix1 is not equal to the rows of MyMatrix2, and can not be multiplied.(MyMatrix::operator*)" << std::endl;
			return res;
		}
		else if (typeid(MyMatrix1.data).name() != typeid(MyMatrix2.data).name()) {
			std::cerr << "Error:Different types of Matrix data(MyMatrix::operator*)" << std::endl;
			return res;
		}
		else
		{
			for (int i = 0; i < MyMatrix1.rows; i++)
				for (int j = 0; j < MyMatrix2.cols; j++)
				{
					ElemType temp = 0.0;
					for (int k = 0; k < MyMatrix1.cols; k++)
					{
						temp += MyMatrix1.data[MyMatrix1.cols * i + k] * MyMatrix2.data[k * MyMatrix2.cols + j];
					}
					res.data[i * res.cols + j] = temp;
				}
			return res;
		}
	}


	//������ + ����(����1)  
	template  <class ElemType>
	MyMatrix<ElemType>  operator+(const MyMatrix<ElemType>& Matrix, const ElemType val)
	{
		MyMatrix<ElemType> res(Matrix.rows, Matrix.cols);
		for (int i = 0; i < Matrix.size; i++)
			res.data[i] = Matrix.data[i] + val;
		return res;
	}
	//������ + ����(����2)
	template  <class ElemType>
	MyMatrix<ElemType>  operator+(const ElemType val, const MyMatrix<ElemType>& Matrix)
	{
		MyMatrix<ElemType> res(Matrix.rows, Matrix.cols);
		for (int i = 0; i < Matrix.size; i++)
			res.data[i] = Matrix.data[i] + val;
		return res;
	}
	//������ - ����(����1)  
	template  <class ElemType>
	MyMatrix<ElemType>  operator-(const MyMatrix<ElemType>& Matrix, const ElemType val)
	{
		MyMatrix<ElemType> res(Matrix.rows, Matrix.cols);
		for (int i = 0; i < Matrix.size; i++)
			res.data[i] = Matrix.data[i] + val;
		return res;
	}
	//������ - ����(����2)
	template  <class ElemType>
	MyMatrix<ElemType>  operator-(const ElemType val, const MyMatrix<ElemType>& Matrix)
	{
		MyMatrix<ElemType> res(Matrix.rows, Matrix.cols);
		for (int i = 0; i < Matrix.size; i++)
			res.data[i] = val - Matrix.data[i];
		return res;
	}
	//������ * ����(����1)  
	template  <class ElemType>
	MyMatrix<ElemType>  operator*(const MyMatrix<ElemType>& Matrix, const ElemType val)
	{
		MyMatrix<ElemType> result(Matrix.rows, Matrix.cols);
		for (int i = 0; i < Matrix.size; i++)
			result.data[i] = Matrix.data[i] * val;
		return result;
	}
	//������ * ����(����2)
	template  <class ElemType>
	MyMatrix<ElemType>  operator*(const ElemType val, const  MyMatrix<ElemType>& Matrix)
	{
		MyMatrix<ElemType> result(Matrix.rows, Matrix.cols);
		for (int i = 0; i < Matrix.size; i++)
			result.data[i] = Matrix.data[i] * val;
		return result;
	}

	//LUP�ֽ�
	template  <class T>
	bool MyMatrix<T>::LUP_Descomposition(MyMatrix<double>& A, MyMatrix<double>& L, MyMatrix<double>& U, MyMatrix<double>& P)
	{
		int N = P.getSize();
		int row = 0;
		for (int i = 0; i<N; i++) {
			P(0, i) = i;
		}
		for (int i = 0; i<N - 1; i++) {
			double p = 0.0;
			for (int j = i; j<N; j++) {
				if (fabs(A(j, i))>p) {
					p = fabs(A(j, i));
					row = j;
				}
			}
			if (0 == p) {
				std::cerr << "Error:MyMatrix singularity and the inverse MyMatrix is incalculable" << std::endl;
				return false;
			}

			//����P[i]��P[row]
			double tmp = P(0, i);
			P(0, i) = P(0, row);
			P(0, row) = tmp;

			double tmp2 = 0;
			for (int j = 0; j<N; j++) {
				//����A[i][j]�� A[row][j]
				tmp2 = A(i, j);
				A(i, j) = A(row, j);
				A(row, j) = tmp2;
			}

			//����ͬLU�ֽ�
			double u = A(i, i), l = 0;
			for (int j = i + 1; j<N; j++) {
				l = A(j, i) / u;
				A(j, i) = l;
				for (int k = i + 1; k<N; k++) {
					A(j, k) = A(j, k) - A(i, k) * l;
				}
			}
		}

		//����L��U
		for (int i = 0; i<N; i++) {
			for (int j = 0; j <= i; j++) {
				if (i != j) {
					L(i, j) = A(i, j);
				}
				else {
					L(i, j) = 1;
				}
			}
			for (int k = i; k<N; k++) {
				U(i, k) = A(i, k);
			}
		}
		return true;
	}
	//LUP������Է���
	template  <class T>
	MyMatrix<double> MyMatrix<T>::LUP_Solve(MyMatrix<double>& L, MyMatrix<double>& U, MyMatrix<double>& P, MyMatrix<double>& b)
	{
		int N = P.getSize();
		MyMatrix<double> x(1, N, 0.0);
		MyMatrix<double> y(1, N, 0.0);
		//�����滻
		for (int i = 0; i < N; i++)
		{
			y(0, i) = b(0, P(0, i));
			for (int j = 0; j < i; j++)
			{
				y(0, i) = y(0, i) - L(i, j) *y(0, j);
			}
		}
		//�����滻
		for (int i = N - 1; i >= 0; i--)
		{
			x(0, i) = y(0, i);
			for (int j = N - 1; j > i; j--)
			{
				x(0, i) = x(0, i) - U(i, j) * x(0, j);
			}
			x(0, i) /= U(i, i);
		}
		return x;
	}
	//LU�ֽ��������
	template  <class T>
	MyMatrix<double> MyMatrix<T>::inverse()
	{
		int N = rows;
		MyMatrix<double> A_mirror(N, N, 0.0);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				A_mirror(i, j) = data[j * rows + i];
		MyMatrix<double> A_inv(N, N, 0.0);
		MyMatrix<double> A_inv_each(1, N, 0.0);
		MyMatrix<double> b(1, N, 0.0);
		MyMatrix<double> L(N, N, 0.0);
		MyMatrix<double> U(N, N, 0.0);
		MyMatrix<double> P(1, N, 0.0);
		if (LUP_Descomposition(A_mirror, L, U, P) == false) {
			A_inv.is_singular = true;
			return A_inv;
		}
		for (int i = 0; i<N; i++)
		{
			//���쵥λ���ÿһ��
			for (int i = 0; i<N; i++)
			{
				b(0, i) = 0;
			}
			b(0, i) = 1;
			A_inv_each = LUP_Solve(L, U, P, b);
			for (int j = 0; j < N; j++)A_inv(j, i) = A_inv_each(0, j);
		}
		return A_inv;
	}
}

#endif 



