#include "commandlineparser.hpp"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>

using namespace std;

CommandLineParser::CommandLineParser ()
	:command(""),
	 parser_string(":h")
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

void CommandLineParser::add_flag_argument(char name, string description) {
	this->arg_to_string[name] = description;
	this->parser_string += name;
	this->flag_to_parameter[name] = false;
}


void CommandLineParser::exactly_one(char name1, char name2) {
	this->exactly_one_arg.push_back(make_pair(name1, name2));
}

void CommandLineParser::not_both(char name1, char name2) {
	this->not_both_arg.push_back(make_pair(name1, name2));
}

void CommandLineParser::parse(int argc, char* argv[]) {
	int c;
	while ((c = getopt(argc, argv, this->parser_string.c_str())) != -1) {
		// check if option exists
		if (this->arg_to_string.find(c) != this->arg_to_string.end()) {
			if (this->flag_to_parameter.find(c) != this->flag_to_parameter.end()) {
				// argument is a flag, set it to true
				this->flag_to_parameter[c] = true;
			} else {
				this->arg_to_parameter[c] = string(optarg);
			}
		} else {
			if (c == 'h') {
				usage();
				throw exception();
			}
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

	// make sure argument combinations are valid
	for (auto it = this->exactly_one_arg.begin(); it != this->exactly_one_arg.end(); ++it) {
		if (this->arg_to_parameter.find(it->first) != this->arg_to_parameter.end()) {
			if (this->arg_to_parameter.find(it->second) != this->arg_to_parameter.end()) {
				ostringstream oss;
				oss << "Error: options -" << it->first << " and -" << it->second << " cannot be used together." << endl;
				throw runtime_error(oss.str());
			}
		} else {
			if (this->arg_to_parameter.find(it->second) == this->arg_to_parameter.end()) {
				ostringstream oss;
				oss << "Error: one of the options -" << it->first << " and -" << it->second << " must be used." << endl;
				throw runtime_error(oss.str());
			}			
		}
	}

	for (auto it = this->not_both_arg.begin(); it != this->not_both_arg.end(); ++it) {
		if ( (this->arg_to_parameter.find(it->first) != this->arg_to_parameter.end()) && (this->arg_to_parameter.find(it->second) != this->arg_to_parameter.end()) ) {
			ostringstream oss;
			oss << "Error: options -" << it->first << " and -" << it->second << " cannot be used together." << endl;
			throw runtime_error(oss.str());			
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

bool CommandLineParser::get_flag(char name) {
	return this->flag_to_parameter.at(name);
}

bool CommandLineParser::exists(char name) {
	return this->arg_to_parameter.find(name) != this->arg_to_parameter.end();
}

void CommandLineParser::usage() {
	cerr << "usage: " << endl << this->command << endl;
	cerr << endl;
	cerr << "options:" << endl;
	for (auto it = this->arg_to_string.begin(); it != this->arg_to_string.end(); ++it) {
		// check if flag
		if (this->flag_to_parameter.find(it->first) != this->flag_to_parameter.end()) {
			cerr << "\t-" << it->first << "\t" << it->second << endl;
			continue;
		} 
		// if there is a default value, get it
		auto d = this->optional.find(it->first);
		if (d != this->optional.end()) {
			if (d->second != ""){
				cerr << "\t-" << it->first << " VAL\t" << it->second << " (default: " << d->second << ")." << endl;
			} else {
				cerr << "\t-" << it->first << " VAL\t" << it->second << "." << endl;
			}
		} else {
			cerr << "\t-" << it->first << " VAL\t" << it->second << "." << endl;
		}
	}
	cerr << endl;
}

void CommandLineParser::info() {
	for (auto it = this->arg_to_string.begin(); it != this->arg_to_string.end(); ++it) {
		if (this->flag_to_parameter.find(it->first) != this->flag_to_parameter.end()) {
			cerr << "-" << (it->first) << "\t" << get_flag(it->first) << endl;
		} else {
			cerr << "-" << (it->first) << "\t" << get_argument(it->first) << endl;
		}
	}
}
