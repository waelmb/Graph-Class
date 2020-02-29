
#include <iostream>
#include "Graph.h"


int main(int argc, char *argv[]){
  graph g;
  std::vector<graph::vertex_label> bfs_rpt;
  std::vector<graph::vertex_label> dfs_rpt;
  bool cycle_found;

  if(argc != 2) 
    std::cout << "usage:  demo2 <filename>\n";
  else {
    if(!g.read_file(argv[1]))
      std::cout << "could not open file '" << argv[1] << "'\n";
  }

  g.dfs(g.name2id("c1"), dfs_rpt, cycle_found);
  std::cout << "\nDFS from vertex '" << "c1" << "' complete\n\n";
  if(cycle_found) 
    std::cout << "   cycle found!\n";
  else
    std::cout << "   no cycle found\n";
  std::cout << "DFS REPORT:\n\n";
  g.disp_report(dfs_rpt);

  if(g.has_cycle()) {
  	std::cout << "g.has_cycle():  true\n";
  }
  else {
  	std::cout << "g.has_cycle():  false\n";
  }
    


  //testing
  
  vector<int> path;
  int dest = g.name2id("c1");
  std::cout << "printing path to "<< g.id2name(dest) << "or " << dest << "\n";	
  bool suc = g.extract_path(dfs_rpt, dest, path);
  
  if(suc) {
  	std::cout << "extract_path() successful\n";
  	
  	for(int i = 0; i < path.size(); i++) {
  		std::cout << path[i] <<" or " << g.id2name(path[i]) << "\n";
	}
	
	printf("end path\n");
  }
  else {
  	printf("extract_path() failed\n");
  }


  return 0;
}

