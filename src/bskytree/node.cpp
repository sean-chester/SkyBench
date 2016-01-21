#include "bskytree/node.h"

void ClearSkyTree(Node& skytree) {
	std::stack<Node> tree_stack;
	PushStack(tree_stack, skytree);

	while (!tree_stack.empty()) {
		tree_stack.top().children.clear();
		tree_stack.pop();
	}
}

void PushStack(std::stack<Node>& tree_stack, Node& skytree) {
	if (skytree.children.size() > 0) {
		tree_stack.push(skytree);

		const uint32_t num_child = skytree.children.size();
		for (unsigned i = 0; i < num_child; i++)
			PushStack(tree_stack, skytree.children[i]);
	}
}
