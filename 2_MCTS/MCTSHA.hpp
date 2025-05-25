#ifndef MCTSHA_HPP
#define MCTSHA_HPP

/* MTCSHA head file */
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <memory>
#include <numeric>
#include <algorithm>
#include <atomic>

using namespace std;
typedef int                 BOOL;
typedef double				DOUBLE;

// class forward
class TreeNode;

template<typename T>
class matrix : public vector<vector<T>> {
public:
    matrix(size_t n, vector<T> m) : vector <vector<T>>(n, vector<T>(m)) {	}
    matrix(size_t n) :vector <vector<T>>(n) {	}
    matrix() {}
    matrix(const matrix<T>& copy) { *this = copy; }
};

template<class T>
vector<T> M_Map(vector<T>& vecTotal, vector<int>& vecId) {
	vector<T> res;
	for (int i = 0; i < vecId.size(); i++) {
		res.push_back(vecTotal.at(vecId.at(i)));
	}
	return res;
}

class Point
{
public:
    Point() {
        x = 0;
        y = 0;
    };

    Point(double _x, double _y) {
        x = _x;
        y = _y;
    };

    Point(const Point& pt)
        : x(pt.x), y(pt.y) {}

    Point& operator = (const Point& pt) {
        this->x = pt.x;
        this->y = pt.y;
        return *this;
    }

    Point operator + (const Point& pt) {
        return Point(this->x + pt.x, this->y + pt.y);
    }

    Point& operator += (const Point& pt) {
        this->x = this->x + pt.x;
        this->y = this->y + pt.y;
        return *this;
    }

    Point operator - (const Point& pt) {
        return Point(this->x - pt.x, this->y - pt.y);
    }

    Point& operator -= (const Point& pt) {
        this->x = this->x - pt.x;
        this->y = this->y - pt.y;
        return *this;
    }

    Point operator * (const int val) {
        return Point(this->x * val, this->y * val);
    }

    Point& operator *= (const int val) {
        this->x *= val;
        this->y *= val;
        return *this;
    }

    Point operator / (const int val) {
        return Point(this->x / val, this->y / val);
    }

    Point& operator /= (const int val) {
        this->x /= val;
        this->y /= val;
        return *this;
    }

    void setPoint(double _x, double _y) {
        x = _x;
        y = _y;
    };

    void setPoint(Point pos) {
        x = pos.x;
        y = pos.y;
    };

    bool operator==(const Point p)
    {
        return x == p.x && y == p.y;
    }

    double x;
    double y;

};

    
// Index reverse search function
template<typename T>
vector<int> M_InvFind(vector<T> vec, T num = 0) {
	vector<int> ret;
	for (size_t i = 0; i < vec.size(); i++) {
		if (vec[i] != num) {
			ret.push_back(i);
		}
	}
	return ret;
}

// Index search function
template<typename T>
vector<int> M_Find(vector<T> vec, T num = 0) {
	vector<int> ret;
	for (size_t i = 0; i < vec.size(); i++) {
		if (vec[i] == num) {
			ret.push_back(i);
		}
	}
	return ret;
}

enum class OptRTn {
	eOptSuccess = 0,					// 优化成功
	eOptFeederParaErr,					// 供料器参数配置错误
	eOptSlotParaErr,					// 供料器槽位不足或参数不一致
	eOptRuntimeOutErr,					// 优化过程超时
	eOptNoMountTaskErr,					// 无贴装任务
	eOptNoHeadAvalErr,					// 无可用贴片头
	eOptTypeConflictErr,				// 元件注册类型和供料器注册类型不一致
	eOptDesigedHeadErr,					// 指定吸嘴头配置无解
	eOptMechLimitErr,					// 受机械限位限制无解
	eOptInterruptErr,					// 用户外部终止优化线程
	eOptCheckParamErr,					// 优化过程参数检查不通过
	eOptUndefinedErr,					// 未定义的错误类型，需要开发人员调试
	eOptExceptionErr,					// 优化过程存在未经处理的异常
	eOptWorkBalanceErr,					// 负载分配无解
	eOptANCAssignmentErr				// 
};

//日志等级
enum class LogLevel {
	SMT_FATAL = 0,
	SMT_ALERT = 100,
	SMT_CRIT = 200,
	SMT_ERROR = 300,
	SMT_WARN = 400,
	SMT_NOTICE = 500,
	SMT_INFO = 600,
	SMT_DEBUG = 700
};

class CommonPara {
public:
	BOOL				HeadInterval = 30;				// Rod Interval
	BOOL				HeadNum;						//  Head Number
    vector<BOOL>		idxHdUsed;						// Useable Head Number
};

class MountPoint {
public:
	double x;
	double y;
	double r;
	int headnum;
};

class OptimizerGroup {
public:
	matrix<double> CpTX;
	matrix<double> CpTY;
	matrix<double> CpTH;
	matrix<int>    CpNo;
	vector<int>	   CpNum;

	void resize(size_t size) {
		CpTX.resize(size), CpTY.resize(size), CpTH.resize(size), CpNo.resize(size);
		CpNum.resize(size, 0);

		for (size_t i = 0; i < size; i++) {
			CpTX[i].clear(), CpTY[i].clear(), CpTH[i].clear(), CpNo[i].clear();
			CpNum[i] = 0;
		}
	}

	void clear() {
		CpTX.clear(), CpTY.clear(), CpTH.clear(), CpNo.clear();
		CpNum.clear();
	}
};

class OptMountPointInfo {
private:
	string			m_strPartName;		// Component name
	string			m_strNozzleName;	// Nozzle name
	Point			m_pPos;				// Mount point position
	DOUBLE			m_dAngle;			// Mount point angle
	DOUBLE			m_dHeight;			// Mount point height
	BOOL			m_BlockID;			// Block ID
	BOOL			m_CircuitID;		// Circuit ID
	BOOL			m_Level;			// Level ID
	BOOL			m_No;				// Mount point number
	BOOL			m_IdentifyNum;		// Mount point identification number
	BOOL			m_Station;			// Station ID
	BOOL			m_Orbit;			// Orbit ID
    BOOL			m_bSkip;			// Skip flag
public:
Point	GetPosition(){
    return m_pPos; 
};
string OptMountPointInfo::GetStrPartName() {
	return m_strPartName;
};
BOOL OptMountPointInfo::GetSkip() {
	return m_bSkip;
};
};

class  OptComponentInfo {
public:
	OptComponentInfo() {
		m_Limit = 0;
		m_HeadOccupy = 0;
		m_Points = 0;
	}
	string			m_strPartName;			// Component name
	string			m_strNozzleName;		// Nozzle name
	string			m_strFeederName;		// Feeder name

	BOOL			m_Limit;
	BOOL			m_HeadOccupy;
	BOOL			m_Points = 0;

	int				GetSlotOccupy(int SlotWidth);
	int				GetLeftSlotOccupy(int SlotWidth);

};

class  MountOptResult {
public:
	matrix<int>		MountCp;		// Placement sequence of all components in the pick-up cycle
	matrix<int>		RodCp;			// Placement sequence of all rods in the pick-up cycle

	void clear();
	void append(int size = 1, int head_num = 6);
};

// ==== RodOptResult ====
class RodOptResult {
public:
	matrix<int> NozzleGroup;				// Rod installation type
	matrix<int> ComponentGroup;				// Rod installation component
	matrix<int> LevelGroup;					// Rod installation level

	vector<int>	Subcycle;					// Subcycle number
	vector<int>	Sequence;					// Subcycle sequence

	void clear();
};

class  Optimizer {
public:
	Optimizer();
    vector<OptimizerGroup> &GetOptimizerGroup();
    vector<OptComponentInfo> &GetComponentInfo();
    CommonPara &GetCommonParaInfo();
    RodOptResult &GetRodOptResult();
    MountOptResult &GetMountOptResult();
    vector<OptMountPointInfo> &GetMountPoint();
 	bool GetOptimizerExit();
   
    double Calculate_Placement_Distance(MountOptResult *pMountOptResult);
	double Calculate_CyclePlacement_Distance(MountOptResult* pMountOptResult, int curCycle);
    void DPMountSequence(vector<MountPoint> &points);
    void Generate_Cycle_Mount_Sequence(RodOptResult *pRodOptResult, MountOptResult *pMountOptResult, int cntCPg, int cntCycle, bool bFastSearch);
    double GetDPSequenceValue(RodOptResult *pRodOptResult, MountOptResult *pMountOptResult, int cntCPg, int cntCycle);

    TreeNode InitRootNode(RodOptResult *pRodOptResult, MountOptResult *pMountOptResult);

    OptRTn TreeSearch_Point_Allocatin(vector<shared_ptr<Optimizer>> m_opt);

    vector<OptimizerGroup>		m_MountPointGroup;			// Component Mount Point Group
    vector<OptComponentInfo>	m_ComponentInfo;			// Component Information
    RodOptResult				m_RodAssignmentResult;		// Rod Assignment Result
    MountOptResult				m_MountAssignmentResult;	// Mount Assignment Result
    CommonPara					m_CommonPara;				// Common Parameters
    vector<OptMountPointInfo>	m_MountPoint;				// Mount Point Information
    std::atomic<BOOL>			m_bThreadExit;				// Thread exit flag
};






#endif