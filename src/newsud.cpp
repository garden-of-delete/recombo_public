#include "newsud.h"

#include <set>

bool newsudStep(char c, int& dx, int& dy, int& dz)
{
   dx = dy = dz = 0;
   switch (c)
   {
      case 'e': dx =  1; return true;
      case 'w': dx = -1; return true;
      case 'n': dy =  1; return true;
      case 's': dy = -1; return true;
      case 'u': dz =  1; return true;
      case 'd': dz = -1; return true;
      default: return false;
   }
}

std::vector<threevector<int> > newsudToVertices(const std::string& s, int x0, int y0, int z0)
{
   std::vector<threevector<int> > vertices;
   int x = x0, y = y0, z = z0;
   vertices.push_back(threevector<int>(x, y, z));
   int dx, dy, dz;
   for (size_t i = 0; i + 1 < s.size(); i++)
   {
      if (!newsudStep(s[i], dx, dy, dz))
         return std::vector<threevector<int> >();
      x += dx;
      y += dy;
      z += dz;
      vertices.push_back(threevector<int>(x, y, z));
   }
   return vertices;
}

bool newsudIsClosed(const std::string& s)
{
   int x = 0, y = 0, z = 0;
   int dx, dy, dz;
   for (size_t i = 0; i < s.size(); i++)
   {
      if (!newsudStep(s[i], dx, dy, dz))
         return false;
      x += dx;
      y += dy;
      z += dz;
   }
   return x == 0 && y == 0 && z == 0;
}

struct ivec3
{
   int x, y, z;
   bool operator<(const ivec3& o) const
   {
      if (x != o.x) return x < o.x;
      if (y != o.y) return y < o.y;
      return z < o.z;
   }
};

bool newsudIsSelfAvoiding(const std::string& s)
{
   std::set<ivec3> visited;
   ivec3 pos = {0, 0, 0};
   visited.insert(pos);
   int dx, dy, dz;
   for (size_t i = 0; i < s.size(); i++)
   {
      if (!newsudStep(s[i], dx, dy, dz))
         return false;
      pos.x += dx;
      pos.y += dy;
      pos.z += dz;
      if (i + 1 < s.size())
      {
         if (!visited.insert(pos).second)
            return false;
      }
   }
   return true;
}
