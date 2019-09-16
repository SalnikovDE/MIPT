#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

class Graph {
public:
    typedef size_t Vertex;
    Graph(size_t _vertex_count, bool _is_directed) {
        is_directed = _is_directed;
        vertex_count = _vertex_count;
        edges_count = 0;
    }

    bool getIsDirected() {
        return is_directed;
    }
    size_t getVertexCount() const {
        return vertex_count;
    }
    size_t getEdgesCount() const {
        return edges_count;
    }

    virtual void add_edge(const Vertex &start, const Vertex &finish) {
        edges_count++;
    }

    virtual std::vector<Vertex> get_all_neighbours(const Vertex &v) const = 0;
protected:
    size_t vertex_count, edges_count;
    bool is_directed;
};

class AdjListGraph: public Graph {
public:
    AdjListGraph(size_t _vertex_count, bool _is_directed):Graph(_vertex_count, _is_directed) {
        Vertexes.resize(vertex_count);
    }
    virtual void add_edge(const Vertex &start, const Vertex &finish) override {
        Graph::add_edge(start, finish);
        Vertexes[start].push_back(finish);
        if(!is_directed) {
            Vertexes[finish].push_back(start);
        }
    }

    size_t GetVertexDeg(Vertex v) {
        return Vertexes[v].size();
    }

    void ReadMatrix (const std::vector<std::vector<int>> &matrix) {
        for (Vertex i = 0; i < matrix.size(); i++){
            for (Vertex j = 0; j < matrix[i].size(); j++) {
                if (matrix[i][j] != 0) {
                    if(j >= i || is_directed) {
                        add_edge(i, j);
                    }
                }
            }
        }
    };
    virtual std::vector<Vertex> get_all_neighbours(const Vertex &v) const override{
        return Vertexes[v];
    }
private:
    std::vector<std::vector<Vertex>> Vertexes;
};
class WeightedAgjListGraph: public AdjListGraph {

};
namespace GraphAlgorithms {
    enum VertexColours{white, gray, black};
    void dfs_find_components(const Graph::Vertex &v, const AdjListGraph &graph,
                             std::vector<VertexColours> &colours, std::vector<Graph::Vertex> &component) {
        colours[v] = gray;
        component.push_back(v);
        for (auto i : graph.get_all_neighbours(v)) {
            if (colours[i] == white) {
                dfs_find_components(i, graph, colours, component);
            }
        }
        colours[v] = black;
    }

    std::vector<std::vector<Graph::Vertex>> find_components(const AdjListGraph & graph) {
        std::vector<std::vector<Graph::Vertex>> result;
        std::vector<VertexColours> colours(graph.getVertexCount(), white);
        for (Graph::Vertex i = 0; i < graph.getVertexCount(); i++) {
            if (colours[i] == white) {
                result.push_back(std::vector<Graph::Vertex> ());
                dfs_find_components(i, graph, colours, result.back());
            }
        }

        return result;
    }


};


int main() {
    int n, k;
    std::vector<std::vector<int>> v;
    std::cin >> n >> k;
    AdjListGraph g (n, false);
    for (int i = 0; i < k; i++) {
        Graph::Vertex start, finish;
        std::cin >> start >> finish;
        g.add_edge(start - 1, finish - 1);
    }
    auto result = GraphAlgorithms::find_components(g);
    std::cout << result.size() << std::endl;

    for (const auto &comp : result) {
        std::cout << comp.size() << std::endl;
        for (const auto &ver : comp) {
            std::cout << ver + 1 << ' ';
        }
        std::cout << std::endl;
    }
}
