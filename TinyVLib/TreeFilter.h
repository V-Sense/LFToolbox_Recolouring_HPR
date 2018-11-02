#ifndef _TREE_FILTER_H
#define _TREE_FILTER_H

#include "Graph.h"

class CTreeFilter
{
public:
	CTreeFilter();
	~CTreeFilter();

	void Init(CMst* tree);
	void SetSigma(float sigma);
	int Filter(float* nodeCosts, int len);
	
private:
	void Clean();
	float* m_weightTbl;
	float* m_tmpCost;
	CMst* m_tree;
};

CTreeFilter::CTreeFilter()
{
	memset(this, 0, sizeof(*this));
}

CTreeFilter::~CTreeFilter()
{
	Clean();
}

void CTreeFilter::Init(CMst* tree)
{
	Clean();
	m_tree = tree;
	m_weightTbl = new float[m_tree->NodeCnt()];
	m_tmpCost = new float[m_tree->NodeCnt()];
}

void CTreeFilter::SetSigma(float sigma)
{
	int numV = m_tree->NodeCnt();
	for (int i = 0; i < numV; i++){
		m_weightTbl[i] = exp(-m_tree->ParentDistance(i) / sigma);
	}
}

int CTreeFilter::Filter(float* nodeCosts, int len)
{
	int numV = m_tree->NodeCnt();
	if (len != numV){
		return -1;
	}
	
	// first step: from leaf to root
	for (int pos = numV - 1; pos >= 0; pos--){ // from last node
		int idx = m_tree->NodeOrder(pos);

		int childNum = m_tree->ChildNum(idx);
		int* childIdx = m_tree->ChildIndex(idx);
		m_tmpCost[idx] = nodeCosts[idx];
		for (int j = 0; j < childNum; j++){
			int idChild = childIdx[j];
			float simi = m_weightTbl[idChild];
			m_tmpCost[idx] += (simi*m_tmpCost[idChild]);
		}
	}

	// second step: from root to leaf
	nodeCosts[0] = m_tmpCost[0];// the first is root node				 
	for (int pos = 1; pos < numV; pos++){ // from second node
		int idx = m_tree->NodeOrder(pos);

		int parent = m_tree->Parent(idx);
		float simi = m_weightTbl[idx];

		nodeCosts[idx] = m_tmpCost[idx] + simi*(nodeCosts[parent] - simi*m_tmpCost[idx]);
	}

	return 0;
}

void CTreeFilter::Clean()
{
	if (m_weightTbl){
		delete[] m_weightTbl;
		m_weightTbl = NULL;
	}
	if (m_tmpCost){
		delete[] m_tmpCost;
		m_tmpCost = NULL;
	}
	m_tree = NULL;
}


#endif