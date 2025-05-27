#include "MCTSHA.hpp"
#include "Tree.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

void ReadGroupFromTxt(const std::string& filename, std::vector<OptimizerGroup>& frontGroup) {
    std::ifstream fgIn(filename);
    std::string line;
    frontGroup = std::vector<OptimizerGroup>(1); // Initialize with one OptimizerGroup
    while (std::getline(fgIn, line)) {
        if (line.find("CpNum:") != std::string::npos) {
            std::istringstream iss(line.substr(line.find(":") + 1));
            int num;
            frontGroup[0].CpNum.clear();
            while (iss >> num) {
                frontGroup[0].CpNum.push_back(num);
            }
        }
        else if (line.find("CpNo:") != std::string::npos) {
            frontGroup[0].CpNo.clear();
            std::string data = line.substr(line.find(":") + 1);
            std::istringstream groupStream(data);
            std::string group;
            while (std::getline(groupStream, group, ';')) {
                std::istringstream iss(group);
                int val;
                std::vector<int> tempVec;
                while (iss >> val) tempVec.push_back(val);
                if (!tempVec.empty()) frontGroup[0].CpNo.push_back(tempVec);
            }
        }
        else if (line.find("CpTH:") != std::string::npos) {
            frontGroup[0].CpTH.clear();
            std::string data = line.substr(line.find(":") + 1);
            std::istringstream groupStream(data);
            std::string group;
            while (std::getline(groupStream, group, ';')) {
                std::istringstream iss(group);
                double val;
                std::vector<double> tempVec;
                while (iss >> val) tempVec.push_back(val);
                if (!tempVec.empty()) frontGroup[0].CpTH.push_back(tempVec);
            }
        }
        else if (line.find("CpTX:") != std::string::npos) {
            frontGroup[0].CpTX.clear();
            std::string data = line.substr(line.find(":") + 1);
            std::istringstream groupStream(data);
            std::string group;
            while (std::getline(groupStream, group, ';')) {
                std::istringstream iss(group);
                double val;
                std::vector<double> tempVec;
                while (iss >> val) tempVec.push_back(val);
                if (!tempVec.empty()) frontGroup[0].CpTX.push_back(tempVec);
            }
        }
        else if (line.find("CpTY:") != std::string::npos) {
            frontGroup[0].CpTY.clear();
            std::string data = line.substr(line.find(":") + 1);
            std::istringstream groupStream(data);
            std::string group;
            while (std::getline(groupStream, group, ';')) {
                std::istringstream iss(group);
                double val;
                std::vector<double> tempVec;
                while (iss >> val) tempVec.push_back(val);
                if (!tempVec.empty()) frontGroup[0].CpTY.push_back(tempVec);
            }
        }
    }
    fgIn.close();
}

void ReadMountPointFromTxt(const std::string& filename, std::vector<OptMountPointInfo>& in_MountPoint) {
    std::ifstream mpIn(filename);
    if (!mpIn.is_open()) {
        std::cerr << "无法打开 MountPoint 文件" << std::endl;
        return;
    }
    in_MountPoint.clear();
    std::string line;
    while (std::getline(mpIn, line)) {
        double x, y;
        string comma = line.substr(0, line.find(","));
        line = line.substr(line.find(",") + 1);
        int no = stoi(line.substr(0, line.find(",")));
        line = line.substr(line.find(",") + 1);
        x = stod(line.substr(0, line.find(",")));
        line = line.substr(line.find(",") + 1);
        y = stod(line.substr(0, line.find(",")));
        line = line.substr(line.find(",") + 1);
        int isSkip = stoi(line.substr(0, line.find(",")));

        OptMountPointInfo pt;
        pt.SetNo(no);
        pt.SetStrPartName(comma);
        pt.SetSkip(isSkip);
        Point pos;
        pos.x = x;
        pos.y = y;
        pt.SetPosition(pos);
        in_MountPoint.push_back(pt);
    }
    mpIn.close();
}

void ReadFrontRodOptResultFromTxt(const std::string& filename, RodOptResult& pRodOptResult) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "无法打开 " << filename << " 文件" << std::endl;
        return;
    }
    pRodOptResult.NozzleGroup.clear();
    pRodOptResult.ComponentGroup.clear();
    pRodOptResult.LevelGroup.clear();
    pRodOptResult.Subcycle.clear();
    pRodOptResult.Sequence.clear();

    std::string line;
    enum Section { NONE, NOZZLE, COMPONENT, LEVEL, SUBCYCLE, SEQUENCE } section = NONE;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line.find("[NozzleGroup]") != std::string::npos) { section = NOZZLE; continue; }
        if (line.find("[ComponentGroup]") != std::string::npos) { section = COMPONENT; continue; }
        if (line.find("[LevelGroup]") != std::string::npos) { section = LEVEL; continue; }
        if (line.find("[Subcycle]") != std::string::npos) { section = SUBCYCLE; continue; }
        if (line.find("[Sequence]") != std::string::npos) { section = SEQUENCE; continue; }

        std::istringstream iss(line);
        if (section == NOZZLE) {
            std::vector<int> row;
            int val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) pRodOptResult.NozzleGroup.push_back(row);
        } else if (section == COMPONENT) {
            std::vector<int> row;
            int val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) pRodOptResult.ComponentGroup.push_back(row);
        } else if (section == LEVEL) {
            std::vector<int> row;
            int val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) pRodOptResult.LevelGroup.push_back(row);
        } else if (section == SUBCYCLE) {
            int val;
            while (iss >> val) pRodOptResult.Subcycle.push_back(val);
        } else if (section == SEQUENCE) {
            int val;
            while (iss >> val) pRodOptResult.Sequence.push_back(val);
        }
    }
    fin.close();
}

int main() {
    Optimizer m_opt;
    ReadGroupFromTxt("Group\\frontGroup.txt", m_opt.frontGroup);
    ReadGroupFromTxt("Group\\rearGroup.txt", m_opt.rearGroup);
    printf("frontGroup size: %zu, rearGroup size: %zu\n", m_opt.frontGroup.size(), m_opt.rearGroup.size());
    vector<pair<int, int>> cpFindIndex;
    std::ifstream fgIn("Group\\cpFindIndex.txt");
    std::string line;
    while (std::getline(fgIn, line)) {
        // 跳过标题行
        if (line.find("Group\\cpFindIndex") != std::string::npos) continue;
        std::istringstream iss(line);
        int first, second;
        char comma;
        if (iss >> first >> comma >> second) {
            cpFindIndex.emplace_back(first, second);
        }
    }
    fgIn.close();
    m_opt.cpFindIndex = cpFindIndex;

    ReadMountPointFromTxt("Group\\MountPoint.txt", m_opt.m_MountPoint);
    ReadMountPointFromTxt("Group\\FrontMountPoint.txt", m_opt.m_front_MountPoint);
    ReadMountPointFromTxt("Group\\RearMountPoint.txt", m_opt.m_rear_MountPoint);
    ReadFrontRodOptResultFromTxt("Group\\FrontRodOptResult.txt", m_opt.m_Front_pRodOptResult);
    ReadFrontRodOptResultFromTxt("Group\\RearRodOptResult.txt", m_opt.m_Rear_pRodOptResult);

    m_opt.TreeSearch_Point_Allocatin();
    system("pause");
    return 0;
}