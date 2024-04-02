#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include "MyMatrix.hpp"
#include <ctime>


using namespace MyMat;
using namespace std;

//对于较大的矩阵,用此函数打印它的缩写
template  <class ElemType>
void print_large_MyMatrix(MyMatrix<ElemType> A)
{
	if ((A.getCols()>1) && (A.getRows()>1)) {
		std::cout << std::setw(10) << A(0, 0) << " , ... , " << std::setw(10) << A(0, A.getCols() - 1) << std::endl;
		std::cout << "    :      , ... ,     :      " << "size=" << A.getRows() << "x" << A.getCols() << std::endl;
		std::cout << std::setw(10) << A(A.getRows() - 1, 0) << " , ... , " << std::setw(10) << A(A.getRows() - 1, A.getCols() - 1) << std::endl;
	}
	else if (A.getRows()>1) {
		std::cout << std::setw(10) << A(0, 0) << std::endl;
		std::cout << "    :      " << "size=" << A.getRows() << "x" << A.getCols() << std::endl;
		std::cout << std::setw(10) << A(A.getRows() - 1, 0) << std::endl;
	}
	else if (A.getCols()>1) {
		std::cout << std::setw(10) << A(0, 0) << " , ... , " << std::setw(10) << A(0, A.getCols() - 1) << "size=" << A.getRows() << "x" << A.getCols() << std::endl;
	}
}


int main(int argc, char *argv[])
{
	cout.setf(ios::showpoint);
	//计时
	clock_t startTime1, endTime1, startTime2, endTime2;
	//double arr[9] = { 1,2,3,4,5,6,7,8,9 };
	//double drr[9] = { 9,8,7,6,5,4,3,2,1 };
	MyMatrix<double> temp(1024, 1024, 3);
	MyMatrix<double> temp1(1024, 1024, 3);
	MyMatrix<double> temp3;
	MyMatrix<double> temp4;
	startTime1 = clock();
	for (int i = 0; i < 5; i++) {
		temp3 = temp * temp1;
	}
	endTime1 = clock();
	cout << "The run time is:" << (double)((endTime1 - startTime1) / 5) << "ms" << endl;
	cout << "-------------------------" << endl;
	double first = double(endTime1 - startTime1);
	startTime2 = clock();
	for (int i = 0; i < 5; i++) {
		temp4 = temp.mul(temp, temp1);
	}
	endTime2 = clock();
	cout << "The time of acceleration:" << (double)((endTime2 - startTime2) / 5) << "ms" << endl;
	double last = double(endTime2 - startTime2);
	cout << "加速比： " << first / last << endl;
	/*cout << temp3 << endl;
	cout << temp4 << endl;*/
	////int *buffer = new int[15];
	////cout << _msize(buffer) << endl;;
	////int a[3];
	////for (int i = 0; i < 3; i++) {
	////	a[i] = i;
	////}

	////MyMatrix<int> temp(2, , 2);
	//
	//int c = 1;
	//MyMatrix<int> temp(2, 2,1);//构造m*n的空矩阵
	//for (int i = 0; i < temp.getRows(); i++) {
	//	for (int j = 0; j < temp.getCols(); j++) {
	//		temp.setValue(i, j, c++);
	//	}
	//}
	////auto M = temp.inverse();
	//int arr[9] = { 1,2,3,4,5,6,7,8,9};
	//MyMatrix<int> temp(arr, 3, 3);
	////temp.setValue(1, 1, 4);
	//MyMatrix<int> temp1(temp);
	//MyMatrix<int> temp2 = temp* temp1;
	//MyMatrix<int> temp3 = temp.mul(temp,temp1);
	//cout << temp<< endl;
	//cout << temp1 << endl;
	//cout << temp2 << endl;
	//cout << temp3 << endl;
	//cout << temp.getValue(1, 0) << endl;
	////MyMatrix<int> temp(5, 3, 1);
	////MyMatrix<int> temp(buffer, m, n);//把buffer转换成m*n矩阵

	//MyMatrix<int> temp2(temp);
	////MyMatrix<int> temp3;
	////temp3 = temp;//所有元素都要拷贝过去，深拷贝

	////temp = temp + temp2;//同样大小的矩阵
	//temp2 = 3 * temp;//每个元素乘以3
	////temp2 = temp - 3;//每个元素加3
	////temp *= temp3;//支持矩阵乘法和矩阵向量乘，如果不能乘返回值
	//MyMatrix<int> temp3 = 3 - temp;
	//cout << temp << endl;//输出temp格式到屏幕
	////cout << temp2 << endl;
	//cout << temp3 << endl;

	//MyMatrix<int> m_a(5, 5, 1);
	//MyMatrix<int> m_b(5, 2);
	//MyMatrix<int> m_c = m_b - m_a;
	//cout << "m_a=" << endl << m_a;
	//cout << "m_b=" << endl << m_b;
	//cout << "m_b-m_a=" << endl << m_c;
	//cout << "m_b+m_a=" << endl << m_b + m_a;
	//cout << "m_b+m_a+m_a=" << endl << m_b + m_a + m_a;
	//cout << "m_b*m_a=" << endl << m_b*m_a;
	//std::default_random_engine e;
	//std::uniform_real_distribution<double> u(0.0, 1.0);
	//MyMatrix<double> m_d(100, 100, 0);
	//MyMatrix<double> m_e(100, 100, 0);
	//for (int i = 0; i<100; i++) {//生成随机矩阵
	//	for (int j = 0; j<100; j++) {
	//		m_d(i, j) = u(e);
	//	}
	//}
	//cout << "m_d=" << endl;
	//print_large_MyMatrix(m_d);
	//m_e = m_d.transpose()*m_d;//构造正定的随机A矩阵
	//cout << "m_e=" << endl;
	//print_large_MyMatrix(m_e);
	//cout << "start to get inverse MyMatrix" << endl;
	//clock_t time_start = clock();
	//MyMatrix<double> m_f = m_e.inverse();
	//std::cout << "time used is " << 1000 * (clock() - time_start) / (double)CLOCKS_PER_SEC << "ms" << std::endl;
	//cout << "m_e.inverse()=" << endl;
	//print_large_MyMatrix(m_f);
	//cout << "m_e*m_e.inverse()=" << endl;
	//print_large_MyMatrix(m_e*m_f);
	system("pause");
	return 0;
}

