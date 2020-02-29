
#include <iostream>
#include "Graph.h"


int main(){
  graph g;
  std::vector<graph::vertex_label> bfs_rpt;

  g._add_edge("a b 1");
  g._add_edge("a c 1");
  g._add_edge("a d ");
  g.add_edge("b", "d");
  //g._add_edge("b d ten");
  //g._add_edge("a");
  /*
  g.add_edge("a", "b");
  g.add_edge("b", "a");

  g.add_edge("b", "c");
  g.add_edge("d", "a");
  g.add_edge("b", "d");
  g.add_edge("a", "b");
  */

  g.display();
  g.bfs(0, bfs_rpt);
  std::cout << "\nBFS REPORT:\n\n";
  g.disp_report(bfs_rpt);
  
  //TESTING PATH
  
  vector<int> path;
  int dest = 2;
  bool suc = g.extract_path(bfs_rpt, dest, path);
  
  if(suc) {
  	printf("extract_path successful.. printing path\n");
  	
  	for(int i = 0; i < path.size(); i++) {
  		printf("%d\n", path[i]);
	}
	
	printf("end path\n");
  }

  return 0;
}

