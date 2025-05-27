#ifndef TREE_HPP
#define TREE_HPP

#include "MCTSHA.hpp"

#include <mutex>
#include <thread>
#include <condition_variable>

class  MountPointInS1 {
public:
	double X = .0;
	double Y = .0;
	double R = .0;
	int No = -1;

};

class NodeState
{
public:
	NodeState();
	double current_value;
	int current_round_index;
	int cumulative_choices;

	bool operator==(const NodeState& other) const;
	double GetCurrentNodeValue();
	void SetCurrentNodeValue(double value);
	int GetCurrentNodeRoundIndex();
	void SetCurrentNodeRoundIndex(int index);
	bool IsTerminal(int round);
	double Calculate_reward();
	NodeState GetNextStateWithChoice(int count);
};

class TreeNode
{
	friend class Optimizer;
private:
	shared_ptr<TreeNode> p_Parent; // parent node
	vector<shared_ptr<TreeNode>> vpChildren_i; // child nodes
	int nVisits_i = 0; // current node visit time
	double totValue_i = 0;
	NodeState state; // current node state
	double valueDis; // placement route for this node

	vector<bool>	m_vecTreeNodePool;
	std::mutex			m_TreeMutex; 
	std::condition_variable			m_Tree; 
	unsigned long hardware_threads; // thread number
public:
	TreeNode();
	TreeNode(shared_ptr<TreeNode> ref);
	TreeNode(const TreeNode& ref);
	void SetNodeState(NodeState in_state);
	NodeState GetNodeState();
	void SetParentNode(shared_ptr<TreeNode> ptr);
	int GetVisitTime();
	void SetVisitTime(int count);
	void VisitTimeAddOne();
	double GetQualityValue();
	void SetQualityValue(double value);
	void QualityValueAdd(double value);
	void AddChildNode(shared_ptr<TreeNode> subNode, shared_ptr<TreeNode> pNode);

	shared_ptr<TreeNode> treePolicy(shared_ptr<TreeNode> node);
	bool DefaultPolicy(shared_ptr<TreeNode> node, double& value);
	shared_ptr<TreeNode> ExpandNode(shared_ptr<TreeNode> node);
	shared_ptr<TreeNode> BestChild(shared_ptr<TreeNode> node, bool isExploration);
	void BackUp(shared_ptr<TreeNode> node, double reward);
	shared_ptr<TreeNode> MCTSFuction(shared_ptr<TreeNode> node);

	bool SingleGenerateMountPointSequence(int count, int idx = 0);

	bool VirtualHIntMethods(int cycle, int dir, int idx);
	bool ScanBaseMethods(int cycle, int dir, int cnt);
	bool NearestMethods(int cycle, int idx);

	int GetRandomSelect();

	MountOptResult				m_mountOptResult;
	RodOptResult				m_RodOptResult;
	vector<OptimizerGroup>		m_mountPointGroup;
	vector<OptMountPointInfo>	m_MountPoint;
	matrix<MountPointInS1>		m_CpT_S1;
	vector<int>					m_cycleList;
	vector<int>					m_childList;
	vector<double>				m_stepSearch;
	vector<double>				m_VirtualRate;
	int							m_leafNodeNum;
	int							m_maxRunRound;
	int							curCycle; // cycle to current node
	int							gantry;	  // gantry to current node 0: front, 1: rear

	int totalCycle;
	double  head_Interval;
	Optimizer* ptr_Optimizer;
	// ** Thread parameters
	std::atomic<int> simulations;
    std::mutex root_mutex;
    std::condition_variable cv;
    std::mutex cv_mutex;
	void MCTSWorker(shared_ptr<TreeNode> node, int num_simulations);
};

#endif