#include <tuple>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "unionFind.h"
#include "getTime.h"

#include "../include/hdbscan/point.h"
#include "../include/hdbscan/armapoint.h"

#include "../include/hdbscan/edge.h"
#include "../include/hdbscan/hdbscan.h"

using namespace std;
using namespace parlay;
using namespace pargeo;

parlay::sequence<pargeo::dendroNode> pargeo::dendrogram(parlay::sequence<pargeo::wghEdge> &edges, size_t n) {
  timer t; t.start();
  std::cout << "d1" << std::endl;
  sequence<pargeo::wghEdge> edgesSorted =
    parlay::sort(make_slice(edges), [&](wghEdge e1, wghEdge e2) {
				      return e1.weight < e2.weight;
				    });
  std::cout << "d2" << std::endl;

  unionFind uf = unionFind<size_t>(n);

  size_t idx = n;

  sequence<size_t> idxMap(n);

  sequence<size_t> sizes(n);

  parallel_for(0, n,[&](int i) {
		      idxMap[i] = i;
		      sizes[i] = 1;
		    });
  std::cout << "d3" << std::endl;

  auto dendro = parlay::sequence<dendroNode>(edgesSorted.size());

  for(size_t i=0; i<n-1; ++i){
    size_t u = uf.find(edgesSorted[i].u);
    size_t v = uf.find(edgesSorted[i].v);
    std::cout << "d31" << std::endl;
    dendro[i] = tuple(idxMap[u], idxMap[v], edgesSorted[i].weight, sizes[u]+sizes[v]);
    uf.link(u, v);
    std::cout << "d32" << std::endl;
    size_t newIdx = uf.find(u);
    std::cout << "d33" << std::endl;
    idxMap[newIdx] = idx;
    sizes[newIdx] = sizes[u]+sizes[v];
    idx++;
  }
  std::cout << "d4" << std::endl;

  cout << "dendrogram-time = " << t.stop() << endl;

  // for (auto d: dendro) {
  //   cout << get<0>(d) << " " << get<1>(d) << " ";
  //   cout << get<2>(d) << " " << get<3>(d) << endl;
  // }
  return dendro;
}
