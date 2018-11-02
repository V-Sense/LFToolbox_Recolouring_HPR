// [8/10/2016 Yinlin.Hu]

#ifndef _HEAP_H_
#define _HEAP_H_

#include <iostream>

class CHeap
{
public:
	CHeap(int size, bool isMinHeap = false); // default is max-heap
	CHeap();

	~CHeap();

	int Init(int size, bool isMinHeap = false);
	int Push(double index, double weight);
	void Push(double* index, double* weight, int count);
	double Pop(double* weight = NULL);
	double Top(double* weight = NULL);
	void Clear();
	int Size();

	void Print();
	void Test();

private:
	template <class T>
	inline void swap(T& a, T& b);
	inline bool prior(double w1, double w2); // is w1 prior to w2?

	bool m_isMinHeap; // min heap or max heap
	double* m_index;
	double* m_weight;
	int m_size;
	int m_validSize;
};

CHeap::CHeap(int size, bool isMinHeap/* = false*/)
{
	memset(this, 0, sizeof(*this));
	Init(size, isMinHeap);
}

CHeap::CHeap()
{
	memset(this, 0, sizeof(*this));
}

CHeap::~CHeap()
{
	if (m_index){
		delete m_index;
	}
	if (m_weight){
		delete m_weight;
	}
}

int CHeap::Init(int size, bool isMinHeap /*= false*/)
{
	m_isMinHeap = isMinHeap;

	if (m_index){
		delete m_index;
	}
	if (m_weight){
		delete m_weight;
	}

	m_index = new double[size];
	m_weight = new double[size];
	m_size = size;
	m_validSize = 0;

	return 0;
}

int CHeap::Push(double index, double weight)
{
	if (m_validSize == m_size){
		printf("PriorityQueue is Full!\n");
		return -1;
	}

	// insert at the last
	m_index[m_validSize] = index;
	m_weight[m_validSize] = weight;
	m_validSize++;

	// adjust the heap from bottom to top
	int i = m_validSize - 1;
	while (prior(m_weight[i], m_weight[(i - 1) / 2])){
		swap(m_weight[i], m_weight[(i - 1) / 2]);
		swap(m_index[i], m_index[(i - 1) / 2]);
		i = (i - 1) / 2; // jump up to the parent
	}

	return i;
}

void CHeap::Push(double* index, double* weight, int count)
{
	for (int i = 0; i < count; i++){
		Push(index[i], weight[i]);
	}
}

double CHeap::Pop(double* weight)
{
	if (m_validSize == 0){
		return -1;
	}

	if (weight){
		*weight = m_weight[0];
	}
	double outIdx = m_index[0];

	// use the last item to overwrite the first
	m_index[0] = m_index[m_validSize - 1];
	m_weight[0] = m_weight[m_validSize - 1];
	m_validSize--;

	// adjust the heap from top to bottom
	double rawIdx = m_index[0];
	double rawWt = m_weight[0];
	int candiPos = 0; // the root
	int i = 1; // left child of the root
	while (i < m_validSize){
		// test right child
		if (i + 1 < m_validSize && prior(m_weight[i + 1], m_weight[i])){
			i++;
		}
		if (prior(rawWt, m_weight[i])){
			break;
		}
		m_index[candiPos] = m_index[i];
		m_weight[candiPos] = m_weight[i];
		candiPos = i;

		i = (i + 1) * 2 - 1; // left child
	}
	m_index[candiPos] = rawIdx;
	m_weight[candiPos] = rawWt;

	return outIdx;
}

double CHeap::Top(double* weight /*= NULL*/)
{
	if (m_validSize == 0){
		return -1;
	}

	if (weight){
		*weight = m_weight[0];
	}
	return m_index[0];
}

void CHeap::Clear()
{
	m_validSize = 0;
}

int CHeap::Size()
{
	return m_validSize;
}

void CHeap::Print()
{
	for (int i = 0; i < m_validSize; i++){
		printf("<%d %.3f> ", m_index[i], m_weight[i]);
	}
	printf("\n");
}

template <class T>
void CHeap::swap(T& a, T& b)
{
	T c = a;
	a = b;
	b = c;
}

bool CHeap::prior(double w1, double w2)
{
	if (m_isMinHeap){
		return w1 < w2;
	}else{
		return w1 > w2;
	}
}

void CHeap::Test()
{
	Push(1, 4.);
	Push(2, 8.);
	Push(3, 3.);
	Push(4, 2.);
	Push(5, 7.);
	Push(6, 6.);

	Print();
	for (int i = 0; i < 6; i++){
		Pop();
		Print();
	}
}

#endif // _HEAP_H_
