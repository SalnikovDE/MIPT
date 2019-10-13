#include <iostream>
#include <vector>


/*
В галактике «Milky Way» на планете «Snowflake» есть N городов,
 некоторые из которых соединены дорогами. Император галактики «Milky Way» решил
 провести инвентаризацию дорог на планете «Snowflake». Но, как оказалось,
 он не силен в математике, поэтому он просит вас сосчитать количество дорог.
 Требуется написать программу, помогающую императору сосчитать количество дорог на планете «Snowflake».
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

    std::vector<Graph::Vertex>
    find_shortest_way_bfs(const AdjListGraph &graph, const Graph::Vertex &start, const Graph::Vertex &finish) {
        std::queue<Graph::Vertex> to_visit;
        std::vector<Graph::Vertex> ancestor(graph.get_vertex_count());
        std::vector<int> distance(graph.get_vertex_count(), -1);
        to_visit.push(start);
        distance[start] = 0;

        while (!to_visit.empty()) {
            auto cur = to_visit.front();
            to_visit.pop();
            for (const auto &n : graph.get_all_neighbours(cur)) {
                if (distance[n] == -1) {
                    distance[n] = distance[cur] + 1;
                    ancestor[n] = cur;
                    to_visit.push(n);
                }
            }
        }
        std::vector<Graph::Vertex> result(distance[finish] + 1);
        auto cur = finish;
        for (int i = result.size() - 1; i >= 0; i--) {
            result[i] = cur;
            cur = ancestor[cur];
        }
        return result;
    }

    std::vector<std::vector<Graph::Vertex>> find_shortest_way_dijkstra(const WeightedAgjListGraph &graph) {}
};

void
add(Graph &g, int i_start, int j_start, int i_finish, int j_finish, const std::vector<std::vector<char>> &world_map,
    int &fake) {
    if (world_map[i_finish][j_finish] == '#') {
        return;
    } else if (world_map[i_finish][j_finish] == '.') {
        g.add_edge(i_start * world_map[0].size() + j_start, i_finish * world_map[0].size() + j_finish);
    } else if (world_map[i_finish][j_finish] == 'W') {
        g.add_edge(i_start * world_map[0].size() + j_start, fake);
        g.add_edge(fake, i_finish * world_map[0].size() + j_finish);
        fake++;
    }
}


int main() {
    int n;
    std::cin >> n;
    AdjListGraph g(n, false);
    std::vector<std::vector<int>> v;
    for (int i = 0; i < n; i++) {
        v.emplace_back(std::vector<int>());
        for (int j = 0; j < n; j++) {
            int a;
            std::cin >> a;
            v[i].push_back(a);
        }
    }
    g.ReadMatrix(v);
    std::cout << g.getEdgesCount();
    return 0;
}
