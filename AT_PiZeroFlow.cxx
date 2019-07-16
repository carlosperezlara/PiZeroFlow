#include <iostream>
#include <TTree.h>
#include "Analysis.h"
#include "AT_PiZeroFlow.h"

AT_PiZeroFlow::AT_PiZeroFlow() : AT_PiZero() {
}

AT_PiZeroFlow::~AT_PiZeroFlow() {
}

void AT_PiZeroFlow::MyInit() {
}

void AT_PiZeroFlow::MyFinish() {
}

void AT_PiZeroFlow::MyExec() {
  // run selector
  AT_PiZero::MyExec();

}
