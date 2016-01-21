#pragma once

#include <stdint.h>

#include <vector>
#include <stack>

#include "common/common.h"

struct Node {
	uint32_t lattice;
	TUPLE point;
	std::vector<Node> children;

	Node(void) {
		lattice = 0;
	}
	Node(uint32_t _lattice): lattice(_lattice) {
	}
};

void ClearSkyTree(Node& skytree);
void PushStack(std::stack<Node>& tree_stack, Node& skytree);
