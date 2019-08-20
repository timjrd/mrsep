#pragma once

#include <climits>
#include <experimental/optional>
#include <utility>
#include <string>
#include <vector>
#include <set>

#include "Base.hh"

class Instance {
public:
  struct StrainType {
    int  id;
    Base base;
  };
  
  struct Position {
    std::vector<StrainType>            strainTypes;
    std::map<int, std::map<int, Base>> readMappings;
  };
  
private:
  std::vector<Position> _positions;
  std::map<int, std::map<int, std::vector<std::pair<int,Base>>>> _readMappings;
  
  Instance() {}

  void updateReadMappings() {
    _readMappings.clear();
    for (int i = 0; i < _positions.size(); i++)
      for (auto const & r : _positions[i].readMappings)
        for (auto const & m : r.second)
          _readMappings[r.first][m.first].push_back({i, m.second});
  }
  
  bool isConsistent() const {
    if (not _positions.empty()) {
      int const strainTypesSize  = _positions.front().strainTypes.size();
      for (Position const & p : _positions)
        if (p.strainTypes.size()  != strainTypesSize)
          return false;
    }
    return true;
  }
  
public:
  // not really const
  std::vector<Position> const & positions() const {
    return _positions;
  }
  
  // not really const
  std::map<int, std::map<int, std::vector<std::pair<int,Base>>>> const & readMappings() const {
    return _readMappings;
  }
  
  bool match(int const strainType, int const read, int const mapping) const {
    for (auto const & p : _readMappings.at(read).at(mapping)) {
      Base const rb = p.second;
      Base const sb = _positions[p.first].strainTypes[strainType].base;
      if (rb != sb)
        return false;
    }
    return true;
  }
  
  Instance filterPositions(std::set<int> const & keep) const {
    Instance result;
    
    for (int const i : keep)
      result._positions.push_back(_positions[i]);
    
    result.updateReadMappings();
    return result;
  }
  
  Instance filterStrainTypes(std::set<int> const & keep) const {
    Instance result;
    
    for (Position const & _this : _positions) {
      result._positions.emplace_back();
      Position & _that = result._positions.back();
      _that.readMappings = _this.readMappings;
      for (int const i : keep)
        _that.strainTypes.push_back(_this.strainTypes[i]);
    }

    result.updateReadMappings();
    return result;
  }
  
  static std::experimental::optional<Instance> fromFile(std::string const & file) {
    std::ifstream stream;
    stream.open(file);
    return fromStream(stream);
  }
  
  static std::experimental::optional<Instance> fromStream(std::istream & stream) {
    Instance result;
    bool failed = false;
    
    try {
      // load strain types
      int id = 0;
      while (not stream.eof()) {
        std::string line;
        std::getline(stream,line);
        
        std::vector<Base> strainType;
        for (char const c : line)
          strainType.push_back(CHAR_TO_BASE.at(c));
        
        if (not strainType.empty()) {
          result._positions.resize(strainType.size());
          
          for (int i = 0; i < strainType.size(); i++)
            result._positions[i].strainTypes.push_back({id, strainType[i]});
          
          id++;
        }
        else
          break;
      }
      
      // load read mappings
      while (not stream.eof()) {
        int read, mapping;
        stream >> read >> mapping;
        
        std::vector<std::pair<Base,int>> readMapping;
        while (not stream.eof() and stream.peek() != '\n') {
          char c;
          int  i;
          stream >> c >> i;
          readMapping.push_back({CHAR_TO_BASE.at(c), i-1});
        }
        
        if (not readMapping.empty())
          for (auto const bi : readMapping)
            result._positions[bi.second].readMappings[read][mapping] = bi.first;
        else
          break;
      }
    }
    catch (std::out_of_range const & e) {
      failed = true;
    }
    
    if (failed or not result.isConsistent())
      return {};
    else {
      result.updateReadMappings();
      return {result};
    }
  }
};
