#include "Tree.hpp"
#include "MCTSHA.hpp"
#include <stack>
#include <random>
#include <cfloat>
#include <mutex>
#include <set>


NodeState::NodeState()
{
    current_value = 0.;
    current_round_index = 0;
    cumulative_choices = 0;
}

bool NodeState::operator==(const NodeState& other) const
{
    return (current_round_index == other.current_round_index) && cumulative_choices == other.cumulative_choices;
}

double NodeState::GetCurrentNodeValue()
{
    return current_value;
}
void NodeState::SetCurrentNodeValue(double value)
{
    current_value = value;
}
int NodeState::GetCurrentNodeRoundIndex()
{
    return current_round_index;
}
void NodeState::SetCurrentNodeRoundIndex(int index)
{
    current_round_index = index;
}
bool NodeState::IsTerminal(int round)
{
    return current_round_index == round;
}

double NodeState::Calculate_reward()
{
    return current_value;
}

NodeState NodeState::GetNextStateWithChoice(int count)
{
    NodeState nextState;
    nextState.SetCurrentNodeValue(0);
    nextState.SetCurrentNodeRoundIndex(current_round_index + 1);
    nextState.cumulative_choices = count;
    return nextState;
}

TreeNode::TreeNode()
{
    p_Parent = nullptr;
    vpChildren_i.clear();
    curCycle = 0;
    m_leafNodeNum = 0;
}

TreeNode::TreeNode(shared_ptr<TreeNode> ref)
{
    this->ptr_Optimizer = ref->ptr_Optimizer;
    this->head_Interval = ref->head_Interval;

    this->gantry = ref->gantry;
    this->m_mountOptResult = ref->m_mountOptResult;
    this->m_mountPointGroup = ref->m_mountPointGroup;
    this->m_RodOptResult = ref->m_RodOptResult;
    this->m_MountPoint = ref->m_MountPoint;
    this->m_CpT_S1 = ref->m_CpT_S1;
    this->m_stepSearch = ref->m_stepSearch;
    this->totalCycle = ref->totalCycle;
    this->curCycle = ref->curCycle;
    this->m_maxRunRound = ref->m_maxRunRound;
    this->m_leafNodeNum = ref->m_leafNodeNum;
    this->m_cycleList = ref->m_cycleList;
    this->m_VirtualRate = ref->m_VirtualRate;
    this->m_childList.resize(m_leafNodeNum);
    std::iota(this->m_childList.begin(), this->m_childList.end(), 0);
}

TreeNode::TreeNode(TreeNode& ref)
{
    this->ptr_Optimizer = ref.ptr_Optimizer;
    this->head_Interval = ref.head_Interval;

    this->gantry = ref.gantry;
    this->m_mountOptResult = ref.m_mountOptResult;
    this->m_mountPointGroup = ref.m_mountPointGroup;
    this->m_RodOptResult = ref.m_RodOptResult;
    this->m_MountPoint = ref.m_MountPoint;
    this->m_CpT_S1 = ref.m_CpT_S1;
    this->m_stepSearch = ref.m_stepSearch;
    this->totalCycle = ref.totalCycle;
    this->curCycle = ref.curCycle;
    this->m_maxRunRound = ref.m_maxRunRound;
    this->m_leafNodeNum = ref.m_leafNodeNum;
    this->m_cycleList = ref.m_cycleList;
    this->m_VirtualRate = ref.m_VirtualRate;
    this->m_childList.resize(m_leafNodeNum);
    std::iota(this->m_childList.begin(), this->m_childList.end(), 0);
}

bool TreeNode::SingleGenerateMountPointSequence(int cycle, int idx)
{
    bool res;
    if (idx  < m_VirtualRate.size())
    {
        res = VirtualHIntMethods(cycle, 0, idx);// 前virRodIntv_num个节点
    }
    else if (idx < m_childList.size() - 1)
    {
        idx -= m_VirtualRate.size();
        int step = ceil((m_stepSearch.at(1) - m_stepSearch.front()) / m_stepSearch.back());
        int dir = idx / step;
        int cntStep = idx % step;
        res = ScanBaseMethods(cycle, dir, cntStep);
    }
    else
    {
        res = NearestMethods(cycle, idx);
    }
    return res;
}

void TreeNode::SetNodeState(NodeState in_state)
{
    state = in_state;
}

NodeState TreeNode::GetNodeState()
{
    return state;
}
void TreeNode::SetParentNode(shared_ptr<TreeNode> ptr)
{
    p_Parent = ptr;
}
int TreeNode::GetVisitTime()
{
    return nVisits_i;
}
void TreeNode::SetVisitTime(int count)
{
    nVisits_i = count;
}
void TreeNode::VisitTimeAddOne()
{
    nVisits_i += 1;
}
double TreeNode::GetQualityValue()
{
    return totValue_i;
}
void TreeNode::SetQualityValue(double value)
{
    totValue_i = value;
}
void TreeNode::QualityValueAdd(double value)
{
    totValue_i = totValue_i + value;
}
void TreeNode::AddChildNode(shared_ptr<TreeNode> subNode, shared_ptr<TreeNode> pNode)
{
    root_mutex.lock();
    subNode->SetParentNode(pNode);
    vpChildren_i.push_back(subNode);
    root_mutex.unlock();
}

bool TreeNode::VirtualHIntMethods(int curCy, int dir, int idx)
{
    int cycle = m_cycleList.at(curCy);
    int idxCPg = 0;
    int curTotalCycle = 0;
    for (size_t cntCpg = 0; cntCpg < m_RodOptResult.Subcycle.size(); cntCpg++) {
        curTotalCycle += m_RodOptResult.Subcycle.at(cntCpg);
        if (curTotalCycle > cycle) {
            idxCPg = cntCpg;
            break;
        }
    }
    double max_dis = DBL_MAX;
    auto resMountPointGroup = m_mountPointGroup;
    auto resMountOptResult = m_mountOptResult;
    
    auto commPara = ptr_Optimizer->GetCommonParaInfo();
   
        vector<int>		tmpSeqCp(commPara.HeadNum, 0);
        vector<double>	tmpOrderR(commPara.HeadNum, -1);
        vector<int>		unsureS;

        matrix<int>     numCpinCpT;

        auto tempMountPointGroup = m_mountPointGroup;
        auto tempMountOptResult = m_mountOptResult;

        for (size_t i = 0; i < tempMountPointGroup.size(); i++) {// mount point group
            numCpinCpT.push_back(tempMountPointGroup.at(i).CpNum);
        }

        for (size_t id = 0; id < m_RodOptResult.ComponentGroup.at(idxCPg).size(); id++)
        {
            if (m_RodOptResult.ComponentGroup.at(idxCPg).at(id) != 0)
            {
                unsureS.push_back(id);
            }
        }
       
        int numUnS = unsureS.size();
        int cntSeqCp = 0;
        double VirtualIntv = m_VirtualRate[idx] * commPara.HeadInterval;
        double lastpointX, lastpointY;
        
        for (int cntUnS = 0; cntUnS < numUnS; cntUnS++) {
            int idxS = unsureS.at(cntUnS);
            int prior = m_RodOptResult.LevelGroup.at(idxCPg).at(idxS);		
            m_CpT_S1.at(idxS).clear();
            for (int cntCpinCpT = 0; cntCpinCpT < numCpinCpT.at(prior).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1); cntCpinCpT++) {
                MountPointInS1 tmpS1;
                tmpS1.X = tempMountPointGroup.at(prior).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1).at(cntCpinCpT) - idxS * VirtualIntv;
                tmpS1.Y = tempMountPointGroup.at(prior).CpTY.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1).at(cntCpinCpT);
                tmpS1.No = tempMountPointGroup.at(prior).CpNo.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1).at(cntCpinCpT);
                m_CpT_S1.at(idxS).push_back(tmpS1);
            }
        }

        while (numUnS) {
            int idxSBest = 0;
            int cntCpinCpTBest = 0;
            int cntUnSBest = 0;
            int priorBest = 0;
            double minD_Cheby = DBL_MAX, minD_Euler = DBL_MAX;
            if (cntSeqCp == 0) {   
                Point leftMountPointPos = { DBL_MAX, DBL_MAX };
                idxSBest = unsureS.front();

                int prior = m_RodOptResult.LevelGroup.at(idxCPg).at(idxSBest);
                int CpinCpTSize = numCpinCpT.at(prior).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1);
                for (int cntCpinCpT = 0; cntCpinCpT < CpinCpTSize; cntCpinCpT++) {
                    double penality_value = 0;
                    int CpIndex = m_CpT_S1.at(idxSBest).at(cntCpinCpT).No - 1;
                    
                    double Pos = m_MountPoint.at(CpIndex).GetPosition().x;
                    if (m_MountPoint.at(CpIndex).GetPosition().x < leftMountPointPos.x) {
                        leftMountPointPos.x = m_MountPoint.at(CpIndex).GetPosition().x + penality_value;
                        leftMountPointPos.y = m_MountPoint.at(CpIndex).GetPosition().y;
                        cntCpinCpTBest = cntCpinCpT;
                    }
                }
                if (abs(leftMountPointPos.y - DBL_MAX) < 1e-3) {
                    unsureS.erase(unsureS.begin() + cntUnSBest);
                    numUnS = numUnS - 1;
                    continue;
                }
                priorBest = prior;
            }
            else {
                for (int cntUnS = 0; cntUnS < numUnS; cntUnS++) {
                    int idxS = unsureS.at(cntUnS);
                    int prior = m_RodOptResult.LevelGroup.at(idxCPg).at(idxS);
                    for (int cntCpinCpT = 0; cntCpinCpT < numCpinCpT.at(prior).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1); cntCpinCpT++) {
                        double D_Cheby = max(abs(m_CpT_S1.at(idxS).at(cntCpinCpT).X - lastpointX), abs(m_CpT_S1.at(idxS).at(cntCpinCpT).Y - lastpointY));
                        double D_Euler = pow((m_CpT_S1.at(idxS).at(cntCpinCpT).X - lastpointX), 2) +
                            pow((m_CpT_S1.at(idxS).at(cntCpinCpT).Y - lastpointY), 2);
                        int CpIndex = m_CpT_S1.at(idxS).at(cntCpinCpT).No - 1;
                        double Pos = m_MountPoint.at(CpIndex).GetPosition().x;
                        
                        if (minD_Cheby > D_Cheby || (abs(minD_Cheby - D_Cheby) < 0.1 && minD_Euler > D_Euler)) {
                            minD_Euler = D_Euler;
                            minD_Cheby = D_Cheby;
                            idxSBest = idxS;
                            cntCpinCpTBest = cntCpinCpT;
                            cntUnSBest = cntUnS;
                            priorBest = prior;
                        }
                    }
                }
                if (abs(minD_Cheby - DBL_MAX) < 1e-3) {
                    string log_message = " ---- function " + string(__FUNCTION__) + " run in Line:" + std::to_string(__LINE__) + " limit ---- ";
                    printf("%s\n", log_message.c_str());
                    return false;
                }
            }

            tempMountOptResult.RodCp.at(cycle).at(idxSBest) = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).No;	// No. of the component in Rodcp
            tmpOrderR.at(idxSBest) = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).R;
            tmpSeqCp.at(commPara.idxHdUsed.at(cntSeqCp)) = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).No;

            // update last point
            lastpointX = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).X;
            lastpointY = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).Y;
            cntSeqCp++;

            // clear allocated point info
            if (tempMountPointGroup.at(priorBest).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).size() != 0) {
                tempMountPointGroup.at(priorBest).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).erase(tempMountPointGroup.at(priorBest).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).begin() + cntCpinCpTBest);
                tempMountPointGroup.at(priorBest).CpTY.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).erase(tempMountPointGroup.at(priorBest).CpTY.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).begin() + cntCpinCpTBest);
                tempMountPointGroup.at(priorBest).CpNo.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).erase(tempMountPointGroup.at(priorBest).CpNo.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).begin() + cntCpinCpTBest);
                tempMountPointGroup.at(priorBest).CpNum.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1)--;
            }
            else {
                string log_message = " ---- function " + string(__FUNCTION__) + " run in Line:" + std::to_string(__LINE__) + " check error ---- ";
                printf("%s\n", log_message.c_str());
                return false;
            }

            vector<int> idxClear = M_Find(m_RodOptResult.ComponentGroup.at(idxCPg), m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest));
            int numClear = idxClear.size();
            for (int cntClear = 0; cntClear < numClear; cntClear++) {
                if (m_RodOptResult.LevelGroup.at(idxCPg).at(idxClear.at(cntClear)) != m_RodOptResult.LevelGroup.at(idxCPg).at(idxSBest)) {
                    continue;
                }
                m_CpT_S1.at(idxClear.at(cntClear)).erase(m_CpT_S1.at(idxClear.at(cntClear)).begin() + cntCpinCpTBest);
            }
            numCpinCpT.at(priorBest).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1) = numCpinCpT.at(priorBest).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1) - 1;
            unsureS.erase(unsureS.begin() + cntUnSBest);
            numUnS = numUnS - 1;
        }
        valueDis = ptr_Optimizer->GetDPSequenceValue(&m_RodOptResult, &tempMountOptResult, idxCPg, cycle);
        if (max_dis > valueDis)
        {
            max_dis = valueDis;
            resMountPointGroup = tempMountPointGroup;
            resMountOptResult = tempMountOptResult;
        }

    if (fabs(valueDis - DBL_MAX) < 1e-3)
    {
        return false;
    }
    valueDis = max_dis;
    root_mutex.lock();
    m_mountPointGroup = resMountPointGroup;
    m_mountOptResult = resMountOptResult;
    root_mutex.unlock();
    return true;
}

bool TreeNode::ScanBaseMethods(int curCy, int dir, int cnt)
{
    int cycle = m_cycleList.at(curCy);
    double min_cycle_distance = DBL_MAX;
    auto commPara = ptr_Optimizer->GetCommonParaInfo();

    int idxCPg = 0;
    int curTotalCycle = 0;
    for (size_t cntCpg = 0; cntCpg < m_RodOptResult.Subcycle.size(); cntCpg++) {
        curTotalCycle += m_RodOptResult.Subcycle.at(cntCpg);
        if (curTotalCycle > cycle) {
            idxCPg = cntCpg;
            break;
        }
    }

    double max_dis = DBL_MAX;
    auto resMountPointGroup = m_mountPointGroup;
    auto resMountOptResult = m_mountOptResult;

    // === Step 3 search direction ===
        vector<int>	head_range ;
        
        for (size_t id = 0; id < m_RodOptResult.ComponentGroup.at(idxCPg).size(); id++)
        {
            if (m_RodOptResult.ComponentGroup.at(idxCPg).at(id) != 0)
            {
                head_range.push_back(id);
            }
        }

        switch (dir) {
        case 0:
            // search from left to right
            break;
        case 1:
            // search from right to left
            reverse(head_range.begin(), head_range.end());
            break;
        case 2:
            // search from the middle
            vector<int> head_range_tmp;
            while (head_range.size()) {
                short middle = head_range.size() / 2;
                head_range_tmp.push_back(head_range.at(middle));
                head_range.erase(head_range.begin() + middle);
            }
            head_range = head_range_tmp;
            break;
        }

        // === Step 4 search begin point === 
            double startPoint = m_stepSearch.front() + cnt * m_stepSearch.back();
            int select_point = -1;

            MountOptResult cur_cycle_result;
            auto tempMountPointGroup = m_mountPointGroup;
            cur_cycle_result = m_mountOptResult;

            vector<Point> assigned_placement_pos(head_range.size(), { 0,0 });
            set<int> tmp_mount_index;

            //  === Step 5 searching ===
            for (short head_index = 0; head_index < head_range.size(); head_index++) {
                short head = head_range.at(head_index);
                int prior = m_RodOptResult.LevelGroup.at(idxCPg).at(head);
                int component_index = m_RodOptResult.ComponentGroup.at(idxCPg).at(head) - 1;
                int bestIndex = -1;
                if (head_index == 0) {
                    double min_horizontal_distance = DBL_MAX;
                    int CpinCpTSize = tempMountPointGroup.at(prior).CpNum.at(component_index);
                    for (int i = 0; i < tempMountPointGroup.at(prior).CpNum.at(component_index); i++) {
                        int mount_index = tempMountPointGroup.at(prior).CpNo.at(component_index).at(i) - 1;
                        if (tmp_mount_index.find(mount_index) != tmp_mount_index.end()) {
                            continue;
                        }
                        double mount_x = m_MountPoint.at(mount_index).GetPosition().x - head * commPara.HeadInterval;
                        double horizontal_distance = abs(mount_x - startPoint);
                        if (min_horizontal_distance > horizontal_distance) {
                            min_horizontal_distance = horizontal_distance;
                            select_point = mount_index;
                            bestIndex = i;
                        }
                    }
                }
                else {
                    select_point = -1;
                    double min_cheby_distance = DBL_MAX;
                    if (cur_cycle_result.RodCp.at(cycle).at(head) != 0) {
                        continue;
                    }
                    for (int i = 0; i < tempMountPointGroup.at(prior).CpNum.at(component_index); i++) {
                        int mount_index = tempMountPointGroup.at(prior).CpNo.at(component_index).at(i) - 1;
                        if (tmp_mount_index.find(mount_index) != tmp_mount_index.end()) {
                            continue;
                        }
                        vector<Point> point_pos;
                        Point pos2firsthead;
                        pos2firsthead.x = m_MountPoint.at(mount_index).GetPosition().x - head * commPara.HeadInterval;
                        pos2firsthead.y = m_MountPoint.at(mount_index).GetPosition().y;

                        
                        point_pos.push_back(pos2firsthead);
                        // fast search
                        for (short next_head_index = 0; next_head_index < head_range.size(); next_head_index++) {
                            short next_head = head_range.at(next_head_index);
                            if (accumulate(cur_cycle_result.RodCp.at(cycle).begin() + next_head, cur_cycle_result.RodCp.at(cycle).end(), 0) == 0) {
                                break;
                            }
                            if (cur_cycle_result.RodCp.at(cycle).at(next_head) == 0) {
                                continue;
                            }
                            point_pos.push_back(assigned_placement_pos.at(next_head_index));
                        }
                        sort(point_pos.begin(), point_pos.end(), [&](Point elem1, Point elem2) {return elem1.x < elem2.x; });
                        double cheby_distance = 0;
                        for (short index = 0; index < point_pos.size() - 1; index++) {
                            cheby_distance += max(abs(point_pos.at(index).x - point_pos.at(index + 1).x), abs(point_pos.at(index).y - point_pos.at(index + 1).y));
                        }
                        if (min_cheby_distance > cheby_distance) {
                            min_cheby_distance = cheby_distance;
                            select_point = mount_index;
                            bestIndex = i;
                        }
                    }
                }

                if (select_point == -1) {
                    continue;
                }
                cur_cycle_result.RodCp.at(cycle).at(head) = select_point + 1;

                assigned_placement_pos.at(head_index).x = m_MountPoint.at(select_point).GetPosition().x - head * commPara.HeadInterval;
                assigned_placement_pos.at(head_index).y = m_MountPoint.at(select_point).GetPosition().y;
                
                tmp_mount_index.insert(select_point);

                // clear the selected mount point
                tempMountPointGroup.at(prior).CpNo.at(component_index).erase(tempMountPointGroup.at(prior).CpNo.at(component_index).begin() + bestIndex);
                tempMountPointGroup.at(prior).CpTX.at(component_index).erase(tempMountPointGroup.at(prior).CpTX.at(component_index).begin() + bestIndex);
                tempMountPointGroup.at(prior).CpTY.at(component_index).erase(tempMountPointGroup.at(prior).CpTY.at(component_index).begin() + bestIndex);
                tempMountPointGroup.at(prior).CpTH.at(component_index).erase(tempMountPointGroup.at(prior).CpTH.at(component_index).begin() + bestIndex);
                tempMountPointGroup.at(prior).CpNum.at(component_index)--;
            }

            if (std::count(m_RodOptResult.ComponentGroup.at(idxCPg).begin(), m_RodOptResult.ComponentGroup.at(idxCPg).end(), 0) != \
                std::count(cur_cycle_result.RodCp.at(cycle).begin(), cur_cycle_result.RodCp.at(cycle).end(), 0)) {
                return false;
            }

            valueDis = ptr_Optimizer->GetDPSequenceValue(&m_RodOptResult, &cur_cycle_result, idxCPg, cycle);
            if (max_dis > valueDis)
            {
                max_dis = valueDis;
                resMountPointGroup = tempMountPointGroup;
                resMountOptResult = cur_cycle_result;
            }
        
            
    if (fabs(valueDis - DBL_MAX) < 1e-3)
    {
        return false;
    }
    valueDis = max_dis;
    root_mutex.lock();
    m_mountPointGroup = resMountPointGroup;
    m_mountOptResult = resMountOptResult;
    root_mutex.unlock();
    return true;
}

bool TreeNode::NearestMethods(int curCy, int idx)
{
    int cycle = m_cycleList.at(curCy);
    int idxCPg = 0;
    int curTotalCycle = 0;
    for (size_t cntCpg = 0; cntCpg < m_RodOptResult.Subcycle.size(); cntCpg++) {
        curTotalCycle += m_RodOptResult.Subcycle.at(cntCpg);
        if (curTotalCycle > cycle) {
            idxCPg = cntCpg;
            break;
        }
    }
    double max_dis = DBL_MAX;
    auto resMountPointGroup = m_mountPointGroup;
    auto resMountOptResult = m_mountOptResult;
    auto commPara = ptr_Optimizer->GetCommonParaInfo();

    vector<int>		tmpSeqCp(commPara.HeadNum, 0);
    vector<double>	tmpOrderR(commPara.HeadNum, -1);
    vector<int>		unsureS;

    matrix<int>     numCpinCpT;

    auto tempMountPointGroup = m_mountPointGroup;
    auto tempMountOptResult = m_mountOptResult;

    for (size_t i = 0; i < tempMountPointGroup.size(); i++) {// calculate the number of different types of components with different priorities
        numCpinCpT.push_back(tempMountPointGroup.at(i).CpNum);
    }

    for (size_t id = 0; id < m_RodOptResult.ComponentGroup.at(idxCPg).size(); id++)
    {
        if (m_RodOptResult.ComponentGroup.at(idxCPg).at(id) != 0)
        {
            unsureS.push_back(id);
        }
    }
    int numUnS = unsureS.size();
    int cntSeqCp = 0;
    double lastpointX, lastpointY;
    // ==== Calculate the distance between the last point and the current point ====
    for (int cntUnS = 0; cntUnS < numUnS; cntUnS++) {
        int idxS = unsureS.at(cntUnS);
        int prior = m_RodOptResult.LevelGroup.at(idxCPg).at(idxS); 
        m_CpT_S1.at(idxS).clear();
        for (int cntCpinCpT = 0; cntCpinCpT < numCpinCpT.at(prior).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1); cntCpinCpT++) {
            MountPointInS1 tmpS1;
            tmpS1.X = tempMountPointGroup.at(prior).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1).at(cntCpinCpT);
            tmpS1.Y = tempMountPointGroup.at(prior).CpTY.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1).at(cntCpinCpT);
            tmpS1.No = tempMountPointGroup.at(prior).CpNo.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1).at(cntCpinCpT);
            m_CpT_S1.at(idxS).push_back(tmpS1);
        }
    }

    // === Rod Task ===
    while (numUnS) {
        int idxSBest = 0;
        int cntCpinCpTBest = 0;
        int cntUnSBest = 0;
        int priorBest = 0;
        double minD_Cheby = DBL_MAX, minD_Euler = DBL_MAX;
        if (cntSeqCp == 0) {// First Mount Point 
            Point leftMountPointPos = { DBL_MAX, DBL_MAX };
            idxSBest = unsureS.front();

            int prior = m_RodOptResult.LevelGroup.at(idxCPg).at(idxSBest);
            for (int cntCpinCpT = 0; cntCpinCpT < numCpinCpT.at(prior).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1); cntCpinCpT++) {
                double penality_value = 0;
                int CpIndex = m_CpT_S1.at(idxSBest).at(cntCpinCpT).No - 1;
                
                if (m_MountPoint.at(CpIndex).GetPosition().x + penality_value < leftMountPointPos.x
                    || (abs(m_MountPoint.at(CpIndex).GetPosition().x - leftMountPointPos.x) < 2 && m_MountPoint.at(CpIndex).GetPosition().y < leftMountPointPos.y)) {
                    leftMountPointPos.x = m_MountPoint.at(CpIndex).GetPosition().x + penality_value;
                    leftMountPointPos.y = m_MountPoint.at(CpIndex).GetPosition().y;
                    cntCpinCpTBest = cntCpinCpT;
                }
            }
            if (abs(leftMountPointPos.y - DBL_MAX) < 1e-3) {
                unsureS.erase(unsureS.begin() + cntUnSBest);
                numUnS = numUnS - 1;
                continue;
            }
            priorBest = prior;
        }
        else {
            for (int cntUnS = 0; cntUnS < numUnS; cntUnS++) {
                int idxS = unsureS.at(cntUnS);
                int prior = m_RodOptResult.LevelGroup.at(idxCPg).at(idxS);
                for (int cntCpinCpT = 0; cntCpinCpT < numCpinCpT.at(prior).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxS) - 1); cntCpinCpT++) {
                    double D_Cheby = max(abs(m_CpT_S1.at(idxS).at(cntCpinCpT).X - lastpointX), abs(m_CpT_S1.at(idxS).at(cntCpinCpT).Y - lastpointY));
                    double D_Euler = pow((m_CpT_S1.at(idxS).at(cntCpinCpT).X - lastpointX), 2) +
                        pow((m_CpT_S1.at(idxS).at(cntCpinCpT).Y - lastpointY), 2);
                    int CpIndex = m_CpT_S1.at(idxS).at(cntCpinCpT).No - 1;
                    double Pos = m_MountPoint.at(CpIndex).GetPosition().x;
                    
                    if (minD_Cheby > D_Cheby || (abs(minD_Cheby - D_Cheby) < 0.1 && minD_Euler > D_Euler)) {
                        minD_Euler = D_Euler;
                        minD_Cheby = D_Cheby;
                        idxSBest = idxS;
                        cntCpinCpTBest = cntCpinCpT;
                        cntUnSBest = cntUnS;
                        priorBest = prior;
                    }
                }
            }
        }

        tempMountOptResult.RodCp.at(cycle).at(idxSBest) = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).No;	// No. in RodCp
        tmpOrderR.at(idxSBest) = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).R;
        tmpSeqCp.at(commPara.idxHdUsed.at(cntSeqCp)) = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).No;

        // update last point
        lastpointX = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).X;
        lastpointY = m_CpT_S1.at(idxSBest).at(cntCpinCpTBest).Y;
        cntSeqCp++;

        // clear allocated component information
        if (tempMountPointGroup.at(priorBest).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).size() != 0) {
            tempMountPointGroup.at(priorBest).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).erase(tempMountPointGroup.at(priorBest).CpTX.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).begin() + cntCpinCpTBest);
            tempMountPointGroup.at(priorBest).CpTY.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).erase(tempMountPointGroup.at(priorBest).CpTY.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).begin() + cntCpinCpTBest);
            tempMountPointGroup.at(priorBest).CpNo.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).erase(tempMountPointGroup.at(priorBest).CpNo.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1).begin() + cntCpinCpTBest);
            tempMountPointGroup.at(priorBest).CpNum.at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1)--;
        }
        else {
            string log_message = " ---- function " + string(__FUNCTION__) + " run in Line:" + std::to_string(__LINE__) + " check error ---- ";
            printf("%s\n", log_message.c_str());
            return false;
        }

        vector<int> idxClear = M_Find(m_RodOptResult.ComponentGroup.at(idxCPg), m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest));
        int numClear = idxClear.size();
        for (int cntClear = 0; cntClear < numClear; cntClear++) {
            if (m_RodOptResult.LevelGroup.at(idxCPg).at(idxClear.at(cntClear)) != m_RodOptResult.LevelGroup.at(idxCPg).at(idxSBest)) {
                continue;	// 
            }
            m_CpT_S1.at(idxClear.at(cntClear)).erase(m_CpT_S1.at(idxClear.at(cntClear)).begin() + cntCpinCpTBest);
        }
        numCpinCpT.at(priorBest).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1) = numCpinCpT.at(priorBest).at(m_RodOptResult.ComponentGroup.at(idxCPg).at(idxSBest) - 1) - 1;
        unsureS.erase(unsureS.begin() + cntUnSBest);
        numUnS = numUnS - 1;
    }
    valueDis = ptr_Optimizer->GetDPSequenceValue(&m_RodOptResult, &tempMountOptResult, idxCPg, cycle);
    if (max_dis > valueDis)
    {
        max_dis = valueDis;
        resMountPointGroup = tempMountPointGroup;
        resMountOptResult = tempMountOptResult;
    }
    if (fabs(valueDis - DBL_MAX) < 1e-3)
    {
        return false;
    }
    valueDis = max_dis;
    root_mutex.lock();
    m_mountPointGroup = resMountPointGroup;
    m_mountOptResult = resMountOptResult;
    root_mutex.unlock();
    return true;
}

int TreeNode::GetRandomSelect()
{
    root_mutex.lock();
    if (m_childList.size() <= 0)
    {
        root_mutex.unlock();
        return m_leafNodeNum;
    }
    std::uniform_int_distribution<int> u(0, m_childList.size() - 1);
    std::mt19937 gen((unsigned int)time(NULL));
    int idx = u(gen);
    int select = m_childList.at(idx);
    m_childList.erase(m_childList.begin() + idx);
    root_mutex.unlock();
    return select;
}

void TreeNode::MCTSWorker(shared_ptr<TreeNode> node, int num_simulations)
{
    for (size_t count = 0; count < num_simulations; count++) {
        if (ptr_Optimizer->GetOptimizerExit())
            return;
        auto expandNode = treePolicy(node);
        double reward = 0;
        bool res = DefaultPolicy(expandNode, reward);
        if (!res)
            return;
        BackUp(expandNode, reward);
    }
    m_TreeMutex.lock();
    m_vecTreeNodePool.push_back(true);
    if (m_vecTreeNodePool.size() == hardware_threads) {
        // Notify the main thread when all threads are done				
        m_Tree.notify_all();
    }
    m_TreeMutex.unlock();
}

bool TreeNode::DefaultPolicy(shared_ptr<TreeNode> node, double& value)
{
    int cycle = node->curCycle;
    shared_ptr<TreeNode> tempNode = make_shared<TreeNode>(node);
    int searchDeep = 5; // simulation depth
    while (1){
        if (ptr_Optimizer->GetOptimizerExit())return false; 
        cycle++;
        if (cycle >= totalCycle || (cycle - node->curCycle) > searchDeep) break;
        int select = node->GetNodeState().cumulative_choices; 
        bool res = tempNode->SingleGenerateMountPointSequence(cycle, select);
        if (!res) break;
        value += tempNode->valueDis;
    }
    return true;
}
shared_ptr<TreeNode> TreeNode::ExpandNode(shared_ptr<TreeNode> node)
{
    shared_ptr<TreeNode> subNode = make_shared<TreeNode>(node);
    subNode->curCycle = node->curCycle + 1;
    int select = node->GetRandomSelect();
    if (select>= m_leafNodeNum)
    {
        return nullptr;
    }
    NodeState newNodeState = node->GetNodeState().GetNextStateWithChoice(select);
    subNode->SingleGenerateMountPointSequence(subNode->curCycle, select);
    newNodeState.current_value = subNode->valueDis;
    subNode->SetNodeState(newNodeState); 
    node->AddChildNode(subNode, node);
    return subNode;
}
void TreeNode::BackUp(shared_ptr<TreeNode> node, double reward)
{
    root_mutex.lock();
    while (node != nullptr)
    {
        node->VisitTimeAddOne();
        node->QualityValueAdd(reward);
        node->totValue_i += node->valueDis;
        node = node->p_Parent;
    }
    root_mutex.unlock();
}
shared_ptr<TreeNode> TreeNode::BestChild(shared_ptr<TreeNode> node, bool isExploration)
{
    double bestScore = DBL_MAX;
    shared_ptr<TreeNode> bestSubNode = nullptr; 
    double parameter = sqrt(2);
    if (node->vpChildren_i.size() <= 0)
    {
        return node;
    }
    root_mutex.lock();
    for (auto subNode : node->vpChildren_i){
        if (isExploration == false)
        {
            parameter = 0; // Select minmun child node after search
            
        }
        if (subNode == nullptr)
        {
            continue;
        }
        double left = subNode->GetQualityValue() / subNode->GetVisitTime();
        double right = log(node->GetVisitTime())/ subNode->GetVisitTime();
        double score = left - parameter * sqrt(right);
        if (score < bestScore){
            bestScore = score;
            bestSubNode = subNode;
        }
    }
    root_mutex.unlock();
    if (bestSubNode == nullptr)
    {
        return node;
    }
    return bestSubNode;
}
shared_ptr<TreeNode> TreeNode::treePolicy(shared_ptr<TreeNode> node)
{
    while (node->curCycle < (node->totalCycle - 1)){
        if (node->m_childList.empty())
            node = BestChild(node, true);
        else {
            auto ptrNode = ExpandNode(node);
            if (ptrNode != nullptr)
            {
                return ptrNode;
            }
        }
    }
    return node;
}

shared_ptr<TreeNode> TreeNode::MCTSFuction(shared_ptr<TreeNode> node)
{
    //hardware_threads = 1;
    hardware_threads = 3;

    std::vector<std::thread> workers;
    int searchCount = 5;
    for (int i = 0; i < hardware_threads; ++i) {
        workers.emplace_back(&TreeNode::MCTSWorker, this, node, searchCount);
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
    // wait for all threads to finish
    // Notify the main thread when all threads are done
    while (m_vecTreeNodePool.size() != hardware_threads) {
        std::unique_lock<mutex> lk(m_TreeMutex);
        m_Tree.wait(lk, [this] {
            return m_vecTreeNodePool.size() == hardware_threads; });
    }

    if (ptr_Optimizer->GetOptimizerExit())
        return nullptr;
    auto child = BestChild(node, false);
    return child;
}