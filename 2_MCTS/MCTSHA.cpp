#include "MCTSHA.hpp"
#include "Tree.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <cassert>

vector<OptimizerGroup>& Optimizer::GetOptimizerGroup()
{
	return m_MountPointGroup;
}

vector<OptComponentInfo>& Optimizer::GetComponentInfo() {
	return m_ComponentInfo;
}

CommonPara& Optimizer::GetCommonParaInfo() {
	return m_CommonPara;
}

RodOptResult& Optimizer::GetRodOptResult() {
	return m_RodAssignmentResult;
}

MountOptResult& Optimizer::GetMountOptResult() {
	return m_MountAssignmentResult;
}

vector<OptMountPointInfo>& Optimizer::GetMountPoint() {
	return m_MountPoint;
}

bool Optimizer::GetOptimizerExit() {
	return m_bThreadExit;
}

double Optimizer::Calculate_Placement_Distance(MountOptResult* pMountOptResult) {
	double totalT_OnPCB = .0;  // Total move distance on PCB
	for (size_t curCycle = 0; curCycle < pMountOptResult->RodCp.size(); curCycle++) {
		vector<double> CompPos_x, CompPos_y, CompPos_r;
		for (size_t h = 0; h < pMountOptResult->MountCp[curCycle].size(); h++) {
			int CpIndex = pMountOptResult->MountCp[curCycle][h];
			int HeadNo;
			if (CpIndex != 0) {
				HeadNo = M_Find(pMountOptResult->RodCp[curCycle], CpIndex).front(); // Find the index of the component in the rod

				CompPos_x.push_back(m_MountPoint[CpIndex - 1].GetPosition().x - HeadNo * m_CommonPara.HeadInterval);
				CompPos_y.push_back(m_MountPoint[CpIndex - 1].GetPosition().y);
			}
		}

		for (size_t i = 1; i < CompPos_x.size(); i++) {
			totalT_OnPCB += max(abs(CompPos_x[i] - CompPos_x[i - 1]), abs(CompPos_y[i] - CompPos_y[i - 1]));
		}
	}
	return totalT_OnPCB;
}

double Optimizer::Calculate_CyclePlacement_Distance(MountOptResult *pMountOptResult, int curCycle)
{
    double totalT_OnPCB = .0;  // Total move distance on PCB
	vector<double> CompPos_x, CompPos_y, CompPos_r;
		for (size_t h = 0; h < pMountOptResult->MountCp[curCycle].size(); h++) {
			int CpIndex = pMountOptResult->MountCp[curCycle][h];
			int HeadNo;
			if (CpIndex != 0) {
				HeadNo = M_Find(pMountOptResult->RodCp[curCycle], CpIndex).front(); // Find the index of the component in the rod

				CompPos_x.push_back(m_MountPoint[CpIndex - 1].GetPosition().x - HeadNo * m_CommonPara.HeadInterval);
				CompPos_y.push_back(m_MountPoint[CpIndex - 1].GetPosition().y);
			}
		}

		for (size_t i = 1; i < CompPos_x.size(); i++) {
			totalT_OnPCB += max(abs(CompPos_x[i] - CompPos_x[i - 1]), abs(CompPos_y[i] - CompPos_y[i - 1]));
		}
	return totalT_OnPCB;
}

void Optimizer::DPMountSequence(vector<MountPoint>& points) {
	vector<MountPoint> sort_points = points;
	MountPoint pick_point;
	pick_point.headnum = -1;
	pick_point.r = 0;
	double points_x_sum = 0;
	for (int i = 0; i < points.size(); i++)
		points_x_sum += points.at(i).x;
	pick_point.x = points_x_sum / points.size();

	pick_point.y = 0;

	sort_points.insert(sort_points.begin(), pick_point);
	vector<vector<vector<int>>> dp_path(1 << (sort_points.size()), matrix<int>(sort_points.size()));
	matrix<double> dp_distance(1 << (sort_points.size()), vector<double>(sort_points.size(), DBL_MAX));
	matrix<double> distance(sort_points.size(), vector<double>(sort_points.size(), DBL_MAX));
    
    
	for (int i = 0; i < sort_points.size() - 1; i++) {
		for (int j = i + 1; j < sort_points.size(); j++) {
			distance[i][j] = distance[j][i] = max(abs(sort_points[i].x - sort_points[j].x), abs(sort_points[i].y - sort_points[j].y));
		}
	}

	int num_pos = sort_points.size();
	for (int s = 1; s < (1 << num_pos); s += 2) {
		// set s include 0
		if (!(s & 1))
			continue;

		for (int j = 1; j < num_pos; j++) {
			// node j must be in s
			if (!(s & (1 << j)))
				continue;

			// set s include 0 and j
            // Initialization
			if (s == ((1 << j) | 1)) {
				dp_path[s][j].push_back(j);
				dp_distance[s][j] = distance[0][j];
			}

			// Update next node
			for (int i = 1; i < num_pos; i++) {
				// node i must not be in s
				if (s & (1 << i))
					continue;

				if (dp_distance[s][j] + distance[j][i] < dp_distance[s | (1 << i)][i]) {
					dp_distance[s | (1 << i)][i] = dp_distance[s][j] + distance[j][i];
					dp_path[s | (1 << i)][i] = dp_path[s][j];
					dp_path[s | (1 << i)][i].push_back(i);
				}
			}
		}
	}
	double min_distance = DBL_MAX;
	vector<int> min_path;
	for (int i = 1; i < num_pos; i++) {
		int s = (1 << num_pos) - 1;
		if (dp_distance[s][i] + distance[i][0] < min_distance) {
			min_path = dp_path[s][i];
			min_distance = dp_distance[s][i] + distance[i][0];
		}
	}
    // Remove all occurrences of 1 from min_path
    for (auto& v : min_path) v -= 1;
	points = M_Map(points, min_path);
}

void Optimizer::Generate_Cycle_Mount_Sequence(RodOptResult* pRodOptResult, MountOptResult* pMountOptResult, int cntCPg, int cntCycle, bool bFastSearch) {
	// === Generate the cycle mount sequence ===
	vector<MountPoint> Points;
	for (size_t idxS = 0; idxS < pMountOptResult->RodCp.at(cntCycle).size(); idxS++) {
		if (pMountOptResult->RodCp.at(cntCycle).at(idxS) == 0) {
			continue;
		}
		MountPoint point;
		int CpIndex = pMountOptResult->RodCp.at(cntCycle).at(idxS) - 1;

		point.x = m_MountPoint.at(CpIndex).GetPosition().x - m_CommonPara.HeadInterval * idxS;
		point.y = m_MountPoint.at(CpIndex).GetPosition().y;

		point.headnum = idxS;
		Points.push_back(point);
	}
	vector<MountPoint> sortPoints;
	// Dynamic Programming
	DPMountSequence(Points);
	sortPoints = Points;
	
	assert(sortPoints.size() == M_InvFind(pRodOptResult->ComponentGroup.at(cntCPg), 0).size());
}

double Optimizer::GetDPSequenceValue(RodOptResult* pRodOptResult, MountOptResult* pMountOptResult, int cntCPg, int cntCycle)
{
	Generate_Cycle_Mount_Sequence(pRodOptResult, pMountOptResult, cntCPg, cntCycle, false);
	double totalDistance = Calculate_CyclePlacement_Distance(pMountOptResult, cntCycle);
	return totalDistance;
}

TreeNode Optimizer::InitRootNode(RodOptResult pRodOptResult)
{
	MountOptResult			tmpMountOptResult;

	int numSubcycle = accumulate(pRodOptResult.Subcycle.begin(), pRodOptResult.Subcycle.end(), 0);
	tmpMountOptResult.MountCp = matrix<int>(numSubcycle, vector<int>(m_CommonPara.HeadNum, 0));
	tmpMountOptResult.RodCp = tmpMountOptResult.MountCp;

	int compGroupsize = pRodOptResult.ComponentGroup.size();
	// ** Generate the component group
	vector<int> cycleList;
	for (int cntCPg = 0; cntCPg < compGroupsize; cntCPg++) {
		int idxCPg = pRodOptResult.Sequence.at(cntCPg) - 1;
		// subcycle range
		int floorCycle, ceilCycle;
		if (idxCPg == 0) {
			floorCycle = 0;
			ceilCycle = pRodOptResult.Subcycle.at(0);
		}
		else {
			floorCycle = accumulate(pRodOptResult.Subcycle.begin(), pRodOptResult.Subcycle.begin() + idxCPg, 0);
			ceilCycle = floorCycle + pRodOptResult.Subcycle.at(idxCPg);
		}
		for (int cntCycle = floorCycle; cntCycle < ceilCycle; cntCycle++) {
			cycleList.push_back(cntCycle);
		}
	}
	TreeNode rootNode;
	rootNode.m_mountOptResult = tmpMountOptResult;
	rootNode.m_RodOptResult = pRodOptResult;
	rootNode.m_CpT_S1 = matrix<MountPointInS1>(m_CommonPara.HeadNum);
	rootNode.totalCycle = numSubcycle;
	rootNode.curCycle = -1;
	rootNode.m_cycleList = cycleList;
	return rootNode;
}

OptRTn Optimizer::TreeSearch_Point_Allocatin()
{
    string log_message = " ---- function run " + string(__FUNCTION__) + " ---- ";
    printf("%s\n", log_message.c_str());
    // 1. vatriables
    int minTotalD_OnPCB = 0;
    MountOptResult			tmpFrontMountOptResult;         // Front MountOptResult
    MountOptResult			tmpRearMountOptResult;          // Rear MountOptResult
    
    vector<OptimizerGroup> tempFrontGroup = frontGroup;
    vector<OptimizerGroup> tempRearGroup = rearGroup;

    int numFrontSubcycle = accumulate(m_Front_pRodOptResult.Subcycle.begin(), m_Front_pRodOptResult.Subcycle.end(), 0);
    int numRearSubcycle = accumulate(m_Rear_pRodOptResult.Subcycle.begin(), m_Rear_pRodOptResult.Subcycle.end(), 0);

    // Use std::make_shared to avoid slicing and ensure proper construction
    std::shared_ptr<TreeNode> frontRootNodePtr = std::make_shared<TreeNode>(InitRootNode(m_Front_pRodOptResult));
    std::shared_ptr<TreeNode> rearRootNodePtr = std::make_shared<TreeNode>(InitRootNode(m_Rear_pRodOptResult));
    TreeNode& frontRootNode = *frontRootNodePtr;
    TreeNode& rearRootNode = *rearRootNodePtr;

    frontRootNode.ptr_Optimizer = this;
    rearRootNode.ptr_Optimizer = this;
    
    double dMaxPosX = max_element(m_MountPoint.begin(), m_MountPoint.end(), [=](OptMountPointInfo& elem1, OptMountPointInfo& elem2)
        {return elem1.GetPosition().x < elem2.GetPosition().x; })->GetPosition().x;
    double dMimPosX = max_element(m_MountPoint.begin(), m_MountPoint.end(), [=](OptMountPointInfo& elem1, OptMountPointInfo& elem2)
        {return elem1.GetPosition().x > elem2.GetPosition().x; })->GetPosition().x;
    double Center = (dMaxPosX + dMimPosX) / 2;
    double Length = dMaxPosX - dMimPosX;
    double LBoundary = Center - Length / 2;
    double RBoundary = Center + Length / 2;

    double SearchStep = (RBoundary - LBoundary) / 6 < 10 ? 10 : (RBoundary - LBoundary) / 6; // SearchStep
    vector<double> m_stepSearch{ LBoundary, RBoundary, SearchStep };
    // Initial VirtualRate
    vector<double> VirtualRate = { 1 / 6.0, 2 / 6.0, 3 / 6.0, 4 / 6.0, 5 / 6.0, 1.0 };

    frontRootNode.m_stepSearch = m_stepSearch;
    rearRootNode.m_stepSearch = m_stepSearch;

    // frontRootNode.m_leafNodeNum = VirtualRate.size() + 3 * ceil((RBoundary - LBoundary) / SearchStep) + 1;
    // rearRootNode.m_leafNodeNum = VirtualRate.size() + 3 * ceil((RBoundary - LBoundary) / SearchStep) + 1;
    frontRootNode.m_leafNodeNum = VirtualRate.size();
    rearRootNode.m_leafNodeNum = VirtualRate.size();
    frontRootNode.m_maxRunRound = 50; // 每次节点搜索深度 roundmin(50, numSubcycle);
    rearRootNode.m_maxRunRound = 50;
    frontRootNode.m_VirtualRate = VirtualRate;
    rearRootNode.m_VirtualRate = VirtualRate;

    frontRootNode.m_mountPointGroup = tempFrontGroup;
    frontRootNode.m_MountPoint = m_MountPoint;
    frontRootNode.gantry = 0;

    rearRootNode.m_mountPointGroup = tempRearGroup;
    rearRootNode.m_MountPoint = m_MountPoint;
    rearRootNode.gantry = 1;

    // Initialize the child node list
    frontRootNode.m_childList.resize(frontRootNode.m_leafNodeNum);
    std::iota(frontRootNode.m_childList.begin(), frontRootNode.m_childList.end(), 0);

    rearRootNode.m_childList.resize(rearRootNode.m_leafNodeNum);
    std::iota(rearRootNode.m_childList.begin(), rearRootNode.m_childList.end(), 0);

    shared_ptr<TreeNode> searchFrontNode = make_shared<TreeNode>(frontRootNode);
    shared_ptr<TreeNode> searchRearNode = make_shared<TreeNode>(rearRootNode);

    std::ofstream outFile;
    outFile.open("rundata(MCTS).txt", std::ios::app);
    outFile << "==== MCTS algothmic run begin ====" << endl;
    clock_t start_time = clock();
    for (size_t cnt = 0; cnt < (numFrontSubcycle + numRearSubcycle); cnt++) {
        printf("==== MCTS algothmic run cnt: %zu ====\n", cnt);
        if (m_bThreadExit)
            return OptRTn::eOptInterruptErr;
        // ** Switch the search gantry
        bool switchFlag = true;
        if ((searchFrontNode->curCycle >= numFrontSubcycle - 1 || searchFrontNode->curCycle > searchRearNode->curCycle) && searchRearNode->curCycle < numRearSubcycle - 1)
        {
            switchFlag = false;
        }

        // Updated MountPoint Info
        shared_ptr<TreeNode> runNode = switchFlag ? make_shared<TreeNode>(searchFrontNode) : make_shared<TreeNode>(searchRearNode);
        runNode = runNode->MCTSFuction(runNode);
        if (runNode != nullptr) {
            auto tempState = runNode->GetNodeState();
            tempState.current_round_index = 0;
            runNode->SetNodeState(tempState);
            // ** update current gantry root node
            int cycle = runNode->m_cycleList.at(runNode->curCycle);
            int idxCPg = 0;
            int curTotalCycle = 0;
            for (size_t cntCpg = 0; cntCpg < runNode->m_RodOptResult.Subcycle.size(); cntCpg++) {
                curTotalCycle += runNode->m_RodOptResult.Subcycle.at(cntCpg);
                if (curTotalCycle > cycle) {
                    idxCPg = cntCpg;
                    break;
                }
            }
            vector<int> curCycleRod2MountResult = runNode->m_mountOptResult.RodCp.at(cycle);
            shared_ptr<TreeNode> updataNode = switchFlag ? make_shared<TreeNode>(searchRearNode) : make_shared<TreeNode>(searchFrontNode);
            for (size_t cntS = 0; cntS < curCycleRod2MountResult.size(); cntS++)
            {
                int cpindex = runNode->m_RodOptResult.ComponentGroup.at(idxCPg).at(cntS) - 1;
                auto cpFind_iter = find_if(cpFindIndex.begin(), cpFindIndex.end(), [&](pair<int, int>a) {return switchFlag ? a.first == cpindex : a.second == cpindex; });
                if (cpFind_iter != cpFindIndex.end())
                {
                    int CpIndex = switchFlag ? cpFind_iter->second : cpFind_iter->first;
                    auto CpNo = updataNode->m_mountPointGroup.at(0).CpNo.at(CpIndex);
                    auto iterErase = find(CpNo.begin(), CpNo.end(), curCycleRod2MountResult[cntS]);
                    if (iterErase == CpNo.end())
                    {
                        continue;
                    }
                    int eraseCnt = iterErase - CpNo.begin();
                    updataNode->m_mountPointGroup.at(0).CpNo.at(CpIndex).erase(updataNode->m_mountPointGroup.at(0).CpNo.at(CpIndex).begin() + eraseCnt);
                    updataNode->m_mountPointGroup.at(0).CpTX.at(CpIndex).erase(updataNode->m_mountPointGroup.at(0).CpTX.at(CpIndex).begin() + eraseCnt);
                    updataNode->m_mountPointGroup.at(0).CpTY.at(CpIndex).erase(updataNode->m_mountPointGroup.at(0).CpTY.at(CpIndex).begin() + eraseCnt);
                    updataNode->m_mountPointGroup.at(0).CpTH.at(CpIndex).erase(updataNode->m_mountPointGroup.at(0).CpTH.at(CpIndex).begin() + eraseCnt);
                    updataNode->m_mountPointGroup.at(0).CpNum.at(CpIndex)--;
                }
            }
            switchFlag ? searchFrontNode = runNode : searchRearNode = runNode;
            switchFlag ? searchRearNode = updataNode : searchFrontNode = updataNode;
        }
        else {
            if (m_bThreadExit)
                return OptRTn::eOptInterruptErr;
            else
                return OptRTn::eOptCheckParamErr;
        }
    }
    double RunTime = clock() - start_time;
    outFile << "==== MCTS algothmic run end ====" << endl;
    outFile << "Timer: " << to_string(RunTime) << endl;

    tmpFrontMountOptResult = searchFrontNode->m_mountOptResult;
    tmpRearMountOptResult = searchRearNode->m_mountOptResult;

    // Recheck the result
    int nTotalPoints_store, nTotalPoints = 0;
    for (size_t cntCycle = 0; cntCycle < tmpFrontMountOptResult.RodCp.size(); cntCycle++) {
        nTotalPoints += M_InvFind(tmpFrontMountOptResult.RodCp.at(cntCycle), 0).size();
    }
    for (size_t cntCycle = 0; cntCycle < tmpRearMountOptResult.RodCp.size(); cntCycle++) {
        nTotalPoints += M_InvFind(tmpRearMountOptResult.RodCp.at(cntCycle), 0).size();
    }
    nTotalPoints_store = nTotalPoints;

    for (auto Point : m_MountPoint) {
        if (Point.GetSkip())
            continue;

        nTotalPoints--;
    }
    if (nTotalPoints != 0) {
        return OptRTn::eOptCheckParamErr;
    }

    // Update the result
    double totalD_OnPCB = 0;
    totalD_OnPCB += Calculate_Placement_Distance(&tmpFrontMountOptResult);
    totalD_OnPCB += Calculate_Placement_Distance(&tmpRearMountOptResult);
    outFile << "MCTS algothmic run value: " << to_string(totalD_OnPCB) << endl;
    outFile.close();

    return OptRTn::eOptSuccess;
}