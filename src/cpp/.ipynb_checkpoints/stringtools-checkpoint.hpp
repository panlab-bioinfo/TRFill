#include <string>
#include <vector>
using namespace std;

std::vector<string> split_find(const string &str,const string &patten) {
	vector<string> res;
	if (str != "") {
		string strs = str + patten;//在字符串末尾也加入分隔符，方便截取最后一段
		size_t pos = strs.find(patten);//find函数返回的是其下标
		while (pos != string::npos) {//使用循环将每一段都存放如vector中
			string temp = strs.substr(0, pos);
			res.push_back(temp);
			strs = strs.substr(pos + 1);//然后拷贝第一个分隔符后面的所有元素，循环
			pos = strs.find(patten);
		}
	}
	return res;
};
