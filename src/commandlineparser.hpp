#ifndef COMMANDLINEPARSER_HPP
#define COMMANDLINEPARSER_HPP

#include <string>
#include <map>
#include <vector>
#include <utility>

class CommandLineParser {
public:
	CommandLineParser();
	void add_command(std::string command);
	void add_mandatory_argument(char name, std::string description);
	void add_optional_argument(char name, std::string default_val, std::string description);
	void add_flag_argument(char name, std::string description);
	/** exactly one of both parameters must be provided **/
	void exactly_one(char name1, char name2);
	/** arguments cannot be used together **/
	void not_both(char name1, char name2);
	void parse(int argc, char* argv[]);
	std::string get_argument(char name);
	bool get_flag(char name);
	bool exists(char name);
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

	std::vector<std::pair<char, char>> exactly_one_arg;
	std::vector<std::pair<char, char>> not_both_arg;
};

#endif // COMMANDLINEPARSER_HPP
