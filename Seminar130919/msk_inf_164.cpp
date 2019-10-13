#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

/*
 Дан неориентированный невзвешенный граф. Для него вам необходимо найти количество вершин,
 лежащих в одной компоненте связности с данной вершиной (считая эту вершину).
 */

class Graph {
public:
    typedef size_t Vertex;

    Graph(size_t _vertex_count, bool _is_directed) {
        is_directed = _is_directed;
        vertex_count = _vertex_count;
        edges_count = 0;
    }

    bool directed() {
        return is_directed;
    }

    virtual size_t get_vertex_deg(const Vertex &v) const = 0;

    size_t get_vertex_count() const {
        return vertex_count;
    }

    size_t get_edges_count() const {
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

class AdjListGraph : public Graph {
public:
    AdjListGraph(size_t _vertex_count, bool _is_directed) : Graph(_vertex_count, _is_directed) {
        Vertexes.resize(vertex_count);
    }

    AdjListGraph(bool _is_directed, const std::vector<std::vector<int>> &matrix) : Graph(
            matrix.size(), _is_directed) {
        Vertexes.resize(vertex_count);
        read_matrix(matrix);
    }

    virtual void add_edge(const Vertex &start, const Vertex &finish) override {
        Graph::add_edge(start, finish);
        Vertexes[start].push_back(finish);
        if (!is_directed) {
            Vertexes[finish].push_back(start);
        }
    }

    virtual size_t get_vertex_deg(const Vertex &v) const override {
        return Vertexes[v].size();
    }


    virtual std::vector<Vertex> get_all_neighbours(const Vertex &v) const override {
        return Vertexes[v];
    }

private:
    std::vector<std::vector<Vertex>> Vertexes;

    void read_matrix(const std::vector<std::vector<int>> &matrix) {
        if (matrix.size() != get_vertex_count()) {
            throw std::out_of_range("wrong matrix size");
        }
        for (Vertex i = 0; i < matrix.size(); i++) {
            for (Vertex j = 0; j < matrix[i].size(); j++) {
                if (matrix[i][j] != 0) {
                    if (j >= i || is_directed) {
                        add_edge(i, j);
                    }
                }
            }
        }
    };
};

class WeightedAgjListGraph : public AdjListGraph {
private:
    std::vector<int> weight;
public:

};
namespace GraphAlgorithms {
    enum VertexColours {
        white, gray, black
    };

    void dfs_find_components(const Graph::Vertex &v, const Graph &graph,
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

    std::vector<std::vector<Graph::Vertex>> find_components(const Graph &graph) {
        std::vector<std::vector<Graph::Vertex>> result;
        std::vector<VertexColours> colours(graph.get_vertex_count(), white);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == white) {
                result.push_back(std::vector<Graph::Vertex>());
                dfs_find_components(i, graph, colours, result.back());
            }
        }

        return result;
    }

    std::vector<std::vector<Graph::Vertex>> find_shortest_way_bfs(const AdjListGraph &graph) {

    }

    std::vector<std::vector<Graph::Vertex>> find_shortest_way_dijkstra(const WeightedAgjListGraph &graph) {}
};



int main() {
    int n, k;
    std::vector<std::vector<int>> v;
    std::cin >> n >> k;
    for (int i = 0; i < n; i++) {
        v.emplace_back(std::vector<int>());
        for (int j = 0; j < n; j++) {
            int state;
            std::cin >> state;
            v[i].push_back(state);
        }
    }

    AdjListGraph g (false, v);
    auto result = GraphAlgorithms::find_components(g);

    for (const auto &i : result) {
        if (std::find(i.begin(), i.end(), k - 1) != i.end()) {
            std::cout << i.size();
            return 0;
        }
    }
}
