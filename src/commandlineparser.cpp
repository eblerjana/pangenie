#include "commandlineparser.hpp"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>

using namespace std;

CommandLineParser::CommandLineParser ()
	:command(""),
	 parser_string(":")
{}

void CommandLineParser::add_command(string command) {
	this->command = command;
}

void CommandLineParser::add_mandatory_argument(char name, string description) {
	this->arg_to_string[name] = description;
	this->parser_string += name;
	this->parser_string += ':';
	this->mandatory.push_back(name);
}

void CommandLineParser::add_optional_argument(char name, string default_val, string description) {
	this->arg_to_string[name] = description;
	this->parser_string += name;
	this->parser_string += ':';
	this->optional[name] = default_val;
}

void CommandLineParser::parse(int argc, char* argv[]) {
	int c;
	while ((c = getopt(argc, argv, this->parser_string.c_str())) != -1) {
		// check if option exists
		if (this->arg_to_string.find(c) != this->arg_to_string.end()) {
			// valid option
			this->arg_to_parameter[c] = string(optarg);
		} else {
			if (c == ':') {
				ostringstream oss;
				oss << "Error: no argument given for option -" << (char) optopt << "."  << endl;
				throw runtime_error(oss.str());
			}
			if (c == '?') {
				ostringstream oss;
				oss << "Error: unknown option -" << (char) optopt << "." << endl;
			}
		}
	}
	// check if mandatory options are all given
	for (auto it = this->mandatory.begin(); it != this->mandatory.end(); ++it) {
		// check if it was seen
		if (this->arg_to_parameter.find(*it) == this->arg_to_parameter.end()) {
			ostringstream oss;
			oss << "Error: option -" << *it << " is mandatory." << endl;
			throw runtime_error(oss.str());
		}
	}
}

string CommandLineParser::get_argument(char name) {
	auto it = this->arg_to_parameter.find(name);
	if (it != this->arg_to_parameter.end()) {
		return this->arg_to_parameter.at(name);
	} else {
		return this->optional[name];
	}
}

void CommandLineParser::usage() {
	cerr << "Usage: " << this->command << endl;
	cerr << endl;
	cerr << "Options:" << endl;
	for (auto it = this->arg_to_string.begin(); it != this->arg_to_string.end(); ++it) {
		cerr << "\t-" << it->first << "=VAL\t" << it->second << endl;
	}
}

void CommandLineParser::info() {
	for (auto it = this->arg_to_string.begin(); it != this->arg_to_string.end(); ++it) {
		cerr << "-" << (it->first) << "\t" << get_argument(it->first) << endl;
	}
}
