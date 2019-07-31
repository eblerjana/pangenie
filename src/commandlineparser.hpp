#ifndef COMMANDLINEPARSER_HPP
#define COMMANDLINEPARSER_HPP

#include <string>
#include <map>
#include <vector>

class CommandLineParser {
public:
	CommandLineParser();
	void add_command(std::string command);
	void add_mandatory_argument(char name, std::string description);
	void add_optional_argument(char name, std::string default_val, std::string description);
	void add_flag_argument(char name, std::string description);
	void parse(int argc, char* argv[]);
	std::string get_argument(char name);
	bool get_flag(char name);
	void usage();
	void info();
private:
	std::string command;
	std::string parser_string;
	std::map<char, std::string> arg_to_string;
	std::map<char, std::string> arg_to_parameter;
	std::map<char, bool> flag_to_parameter;
	std::vector<char> mandatory;
	std::map<char,std::string> optional;
};

#endif // COMMANDLINEPARSER_HPP
