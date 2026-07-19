#ifndef STRIP_MAP_H
#define STRIP_MAP_H

#include <array>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

static constexpr int kNX = 16;
static constexpr int kNY = 16;
static constexpr int kMaxBoard = 8;
static constexpr int kMaxChannel = 64;

enum class StripAxis { kNone = 0, kX = 1, kY = 2 };

struct MappedStrip {
  StripAxis axis = StripAxis::kNone;
  int strip = -1;
  bool valid = false;
};

struct StripMap {
  std::array<std::array<MappedStrip, kMaxChannel>, kMaxBoard> table{};

  bool Lookup(int board, int channel, MappedStrip &out) const {
    if (board < 0 || board >= kMaxBoard || channel < 0 || channel >= kMaxChannel) {
      return false;
    }
    out = table[board][channel];
    return out.valid;
  }

  bool IsX(int board, int channel) const {
    MappedStrip m;
    return Lookup(board, channel, m) && m.axis == StripAxis::kX;
  }

  bool IsY(int board, int channel) const {
    MappedStrip m;
    return Lookup(board, channel, m) && m.axis == StripAxis::kY;
  }
};

static StripAxis ParseAxisToken(const std::string &token) {
  if (token.empty()) {
    return StripAxis::kNone;
  }
  const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
  if (c == 'X') {
    return StripAxis::kX;
  }
  if (c == 'Y') {
    return StripAxis::kY;
  }
  return StripAxis::kNone;
}

static bool LoadStripMap(const char *path, StripMap &map) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Cannot open channel map: " << path << std::endl;
    return false;
  }

  for (auto &row : map.table) {
    for (auto &cell : row) {
      cell = MappedStrip{};
    }
  }

  std::string line;
  int lineNumber = 0;
  int loaded = 0;
  while (std::getline(in, line)) {
    ++lineNumber;
    if (line.empty()) {
      continue;
    }
    const auto hash = line.find('#');
    if (hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    if (line.find_first_not_of(" \t") == std::string::npos) {
      continue;
    }

    std::istringstream iss(line);
    int board = -1;
    int channel = -1;
    std::string axisToken;
    int strip = -1;
    if (!(iss >> board >> channel >> axisToken >> strip)) {
      std::cerr << path << ":" << lineNumber << ": expected: board channel axis strip" << std::endl;
      return false;
    }

    const StripAxis axis = ParseAxisToken(axisToken);
    if (axis == StripAxis::kNone) {
      std::cerr << path << ":" << lineNumber << ": axis must be X or Y" << std::endl;
      return false;
    }
    if (board < 0 || board >= kMaxBoard || channel < 0 || channel >= kMaxChannel) {
      std::cerr << path << ":" << lineNumber << ": board/channel out of range" << std::endl;
      return false;
    }
    if (axis == StripAxis::kX && (strip < 0 || strip >= kNX)) {
      std::cerr << path << ":" << lineNumber << ": X strip out of range" << std::endl;
      return false;
    }
    if (axis == StripAxis::kY && (strip < 0 || strip >= kNY)) {
      std::cerr << path << ":" << lineNumber << ": Y strip out of range" << std::endl;
      return false;
    }

    map.table[board][channel].axis = axis;
    map.table[board][channel].strip = strip;
    map.table[board][channel].valid = true;
    ++loaded;
  }

  std::cout << "Loaded " << loaded << " channel-map entries from " << path << std::endl;
  return loaded > 0;
}

#endif  // STRIP_MAP_H
