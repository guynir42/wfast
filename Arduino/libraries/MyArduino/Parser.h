#ifndef PARSER_H
#define PARSER_H

#if defined(__AVR_ATmega1280__) || defined(__AVR_ATmega2560__)
#define NUM_COMMANDS 20
#else
#define NUM_COMMANDS 6
#endif

#include <Arduino.h>

typedef void (*FunctionPointer)(char*);

typedef struct {
 
  char keyword[16];
  FunctionPointer fp;
  
} Command;

class Parser {
  
public:
  
  Parser();
  Parser(char *name);
  Parser(char *name, char *version);
  ~Parser();
  
  void printout();
  static void print(int number);
  static void print(char *str);
  static void print(const char *str);
  static void println(int number);
  static void println(char *str);
  static void println(const char *str);
  
  // setup
  void addCommand(char *keyword, FunctionPointer fp);
  
  // getters
  bool isStringComplete() const;
  char inputString() const;
  
  void getVersion() const;
  void getThisSensorName() const;
  
  // setters 
  void setVersion(char *name);
  void setThisSensorName(char *name);
  
  // parsing text
  void readingInput();
  void checkInput();
  void parse(char *text);
  
  static void splitStr(char *input_str, char *out_first_word);
  static void splitStr(char *input_str, char *out_first_word, char *out_second_word);
  
  static int partialMatch(char *str, char *comp_str);
  
  static void itoa(char *buf, int number);
  static void ftoa(char *buf, float number, int precision=4);
    
  static String itos(int number);
  static String ftos(float number, int precision=4);
  
  static int parseBool(char *str);
  
  static void lowerString(char *str);
    
  static void cleanString(char *str);
  
//   static String separators;
  static char separators[3];
  static char line_ending;
#if defined(__AVR_ATmega1280__) || defined(__AVR_ATmega2560__)
  static const int STRN=32;  
  static const int MAXN=128;
#else
  static const int STRN=20;  
  static const int MAXN=36;
#endif  
  static bool debug_bit;
  char sensor_name[STRN];
  
  
protected:
  
  //int _default_num_commands;
  
  void initialize();
  void allocateSpace();
  void clearCommands();
  
  int _num_commands;
  Command _commands[NUM_COMMANDS];
  //int _mem_size;
  
  char _input_string[MAXN];
  int  _input_pos;
  bool _is_done_reading;
  char _input_char;
  
  char _version[MAXN];
  
};

#endif
