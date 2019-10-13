#include <iostream>
#include <vector>
/*
В подземелье M тоннелей и N перекрестков, каждый тоннель соединяет какие-то два перекрестка.
 Мышиный король решил поставить по светофору в каждом тоннеле перед каждым перекрестком. Напишите программу,
 которая посчитает, сколько светофоров должно быть установлено на каждом из перекрестков.
 Перекрестки пронумерованы числами от 1 до N.
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
    int n, m;
    std::cin >> n >> m;
    AdjListGraph v(n, false);
    for (int i = 0; i < m; i++) {
        Graph::Vertex s, f;
        std::cin >> s >> f;
        v.addEdge(s - 1, f - 1);
    }

    for (int i = 0; i < n; i++) {
        std::cout << v.GetVertexDeg(i) << std::endl;
    }
    return 0;
}

