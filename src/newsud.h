#ifndef NEWSUD_H
#define NEWSUD_H

#include <string>
#include <vector>

#include "threevector.h"

bool newsudStep(char c, int& dx, int& dy, int& dz);
std::vector<threevector<int> > newsudToVertices(const std::string& s, int x0 = 0, int y0 = 0, int z0 = 0);

bool newsudIsClosed(const std::string& s);
bool newsudIsSelfAvoiding(const std::string& s);

#endif /* NEWSUD_H */
