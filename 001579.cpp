#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <exception>
#include <cstdlib>
#include <string>
#include <set>
#include <map>
#include <unordered_map>

/*
 * Дан ориентированный граф, в котором могут быть кратные ребра и петли. Каждое ребро имеет вес, выражающийся целым числом (возможно, отрицательным). Гарантируется, что циклы отрицательного веса отсутствуют.

Требуется посчитать длины кратчайших путей от вершины номер 1 до всех остальных вершин.
*/

template<class T, class Hash = std::hash<T>>
class DisjointSet {
public:
    class DSUException {
    };

    DisjointSet(std::vector<T> _v) {
        for (const T &elem : _v) {
            size[elem] = 1;
            prev[elem] = elem;
        }
    }

    DisjointSet(T start, T finish) { // for iterable type, not a reference since copy is needed for cycle
        for (; start <= finish; start++) {
            size[start] = 1;
            prev[start] = start;
        }
    }

    virtual void union_sets(const T &first, const T &second) {
        T first_root = get_set_id(first);
        T second_root = get_set_id(second);
        if (first_root == second_root) {
            throw DSUException();
        }
        if (size[second_root] > size[first_root]) {
            // in case std::swap is undefined for T
            T tmp = first_root;
            first_root = second_root;
            second_root = tmp;
        }
        prev[second_root] = first_root;
        size[first_root] += size[second_root];

    }

    size_t get_set_size(const T &elem) {
        return size[get_set_id(elem)];
    }

    bool in_same_set(const T &first, const T &second) {
        return get_set_id(first) == get_set_id(second);
    }


protected:
    std::unordered_map<T, size_t, Hash> size;
    std::unordered_map<T, T, Hash> prev;

    T get_set_id(const T &elem) {
        if (elem == prev[elem]) {
            return elem;
        }
        return prev[elem] = get_set_id(prev[elem]);
    }
};

template<class T, class Hash = std::hash<T>>
class WeightedDisjointSet : public DisjointSet<T, Hash> {
public:
    WeightedDisjointSet(T start, T finish) : DisjointSet<T, Hash>(start, finish) {

    }

    virtual void union_sets(const T &first, const T &second, int w = 0) {
        auto res_weight = weights[this->get_set_id(first)] + weights[this->get_set_id(second)] + w;
        weights[this->get_set_id(first)] = res_weight;
        weights[this->get_set_id(second)] = res_weight;
        DisjointSet<T, Hash>::union_sets(first, second);
    }

    void add_weight(const T &elem, int w) {
        weights[this->get_set_id(elem)] += w;
    }

    int get_weight(const T &elem) {
        return weights[this->get_set_id(elem)];
    }

private:
    std::unordered_map<T, long long int, Hash> weights;
};

class Graph {
public:
    typedef size_t Vertex;

    Graph(size_t _vertex_count, bool _is_directed) {
        is_directed = _is_directed;
        vertex_count = _vertex_count;
        edges_count = 0;
    }

    bool directed() const {
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

    void add_edge(const Vertex &start, const Vertex &finish) override {
        Graph::add_edge(start, finish);
        Vertexes[start].push_back(finish);
        if (!is_directed) {
            Vertexes[finish].push_back(start);
        }
    }

    size_t get_vertex_deg(const Vertex &v) const override {
        return Vertexes[v].size();
    }


    std::vector<Vertex> get_all_neighbours(const Vertex &v) const override {
        return Vertexes[v];
    }

private:
    std::vector<std::vector<Vertex>> Vertexes;

    void read_matrix(const std::vector<std::vector<int>> &matrix) {
        if (matrix.size() != get_vertex_count()) {
            throw std::out_of_range("wrong matrix size");
        }
        for (Vertex i = 0; i < matrix.size(); i++) {
            if (matrix[i].size() != get_vertex_count()) {
                throw std::out_of_range("wrong matrix size");
            }
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

class AdjMatrixGraph : public Graph {
public:
    AdjMatrixGraph(size_t _vertex_count, bool _is_directed) : Graph(_vertex_count, _is_directed) {
        matrix.resize(vertex_count, std::vector<int>(vertex_count, 0));
    }
    void add_edge(const Vertex &start, const Vertex &finish) override {
        Graph::add_edge(start, finish);
        matrix[start][finish] = 1;
        if (!is_directed) {
            matrix[finish][start] = 1;
        }
    }

    size_t get_vertex_deg(const Vertex &v) const override {
        size_t sum = 0;
        for (Vertex i = 0; i < vertex_count; i++) {
            sum += matrix[v][i];
        }
    }


    std::vector<Vertex> get_all_neighbours(const Vertex &v) const override {
        std::vector<Vertex> neighbours;
        for (Vertex i = 0; i < vertex_count; i++) {
            if (matrix[v][i]) {
                neighbours.push_back(i);
            }
        }
        return neighbours;
    }


protected:
    std::vector<std::vector<int>> matrix;
};

class WeightedAgjListGraph : public AdjListGraph {
public:
    WeightedAgjListGraph(size_t _vertex_count, bool _is_directed) : AdjListGraph(_vertex_count, _is_directed) {

    }

    void add_weighted_edge(const Graph::Vertex &start, const Graph::Vertex &finish, int w)  {
        add_edge(start, finish);
        if (weights.find({start, finish}) != weights.end()) {
            weights[{start, finish}] = std::min(weights[{start, finish}], w);
        }
        weights[{start, finish}] = w;
        if (!is_directed) {
            weights[{finish, start}] = std::min(weights[{start, finish}], weights[{finish, start}]);
        }
    }

    int get_weight(const Graph::Vertex &start, const Graph::Vertex &finish) const  {
        if (weights.find({start, finish}) == weights.end()) {
            return INT16_MAX;
        }
        return weights.at({start, finish});
    }

    virtual std::map<std::pair<Graph::Vertex, Graph::Vertex>, int> get_all_edges() const  {
        return weights;
    }

private:
    std::map<std::pair<Graph::Vertex, Graph::Vertex>, int> weights;
};

class WeightedAdjMatrixGraph : public AdjMatrixGraph {
public:
    WeightedAdjMatrixGraph (size_t _vertex_count, bool _is_directed) : AdjMatrixGraph(_vertex_count, _is_directed) {
        weights.resize(vertex_count, std::vector<int>(vertex_count, INT16_MAX));
        for (int i = 0; i < vertex_count; i++) {
            weights[i][i] = 0;
        }
    }

    void add_weighted_edge(const Graph::Vertex &start, const Graph::Vertex &finish, int w)  {
        add_edge(start, finish);
        weights[start][finish]  = std::min (w, weights[start][finish]);
        if (!is_directed) {
            weights[finish][start] = weights[start][finish];
        }
    }

    int get_weight(const Graph::Vertex &start, const Graph::Vertex &finish) const  {
        return weights[start][finish];
    }
private:
    std::vector<std::vector<int>> weights;
};



namespace GraphAlgorithms {
    enum vertexColours {
        WHITE, GRAY, BLACK
    };
    typedef std::vector<Graph::Vertex> Component;
    typedef std::pair<Graph::Vertex, Graph::Vertex> Edge;

    AdjListGraph get_transposed_adj_list_graph(const Graph &graph) {
        AdjListGraph result_graph(graph.get_vertex_count(), graph.directed());
        for (Graph::Vertex start = 0; start < graph.get_vertex_count(); start++) {
            for (const Graph::Vertex &finish : graph.get_all_neighbours(start)) {
                result_graph.add_edge(finish, start);
            }
        }
        return result_graph;
    }

    void dfs_graph_top_sort(Graph::Vertex ver, const Graph &graph, std::vector<vertexColours> &colours,
                            std::vector<Graph::Vertex> &top_sort_result) {
        colours[ver] = GRAY;
        for (const Graph::Vertex &v : graph.get_all_neighbours(ver)) {
            if (colours[v] == WHITE) {
                dfs_graph_top_sort(v, graph, colours, top_sort_result);
            }
        }
        colours[ver] = BLACK;
        top_sort_result.push_back(ver);

    }

    bool dfs_bipartite(const Graph::Vertex &v, const Graph &graph, std::vector<vertexColours> &colours,
                       std::vector<int> &part) {
        colours[v] = GRAY;
        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == GRAY) {
                if (part[i] == part[v]) {
                    return false;
                }
            } else if (colours[i] == WHITE) {
                part[i] = -part[v];
                if (!dfs_bipartite(i, graph, colours, part)) {
                    return false;
                }
            }
        }
        colours[v] = BLACK;
        return true;
    }

    void dfs_find_components(const Graph::Vertex &v, const Graph &graph,
                             std::vector<vertexColours> &colours, Component &component) {
        colours[v] = GRAY;
        component.push_back(v);
        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == WHITE) {
                dfs_find_components(i, graph, colours, component);
            }
        }
        colours[v] = BLACK;
    }

    bool dfs_check_cicle(const Graph::Vertex &v, const Graph &graph,
                         std::vector<vertexColours> &colours) {
        colours[v] = GRAY;
        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == WHITE) {
                if (dfs_check_cicle(i, graph, colours)) return true;
            } else if (colours[i] == GRAY) {
                colours[v] = BLACK;
                return true;
            }
        }
        colours[v] = BLACK;
        return false;
    }


    Graph::Vertex dfs_find_cicle(const Graph::Vertex &v, const Graph &graph,
                                 std::vector<vertexColours> &colours, std::vector<Graph::Vertex> &ancestors) {
        colours[v] = GRAY;
        Graph::Vertex last_vertex;
        const Graph::Vertex NOT_FOUND = graph.get_vertex_count();

        for (const Graph::Vertex &i : graph.get_all_neighbours(v)) {
            if (colours[i] == GRAY) {
                colours[v] = BLACK;
                ancestors[i] = v;
                return i;
            } else if (colours[i] == WHITE) {
                ancestors[i] = v;
                last_vertex = dfs_find_cicle(i, graph, colours, ancestors);
                if (last_vertex != NOT_FOUND) return last_vertex;
            }
        }
        colours[v] = BLACK;
        return NOT_FOUND;
    }

    void dfs_find_briges(const Graph::Vertex &v, const Graph::Vertex &ancestor, const Graph &graph,
                         std::vector<vertexColours> &colours, std::set<Edge> &bridges, size_t &timer,
                         std::vector<int> &tin, std::vector<int> &tup) {
        colours[v] = GRAY;
        tup[v] = tin[v] = timer++;
        for (const Graph::Vertex &to : graph.get_all_neighbours(v)) {
            if (to == ancestor) {
                continue;
            } else if (colours[to] == GRAY) {
                tup[v] = std::min(tup[v], tin[to]);
            } else if (colours[to] == WHITE) {
                dfs_find_briges(to, v, graph, colours, bridges, timer, tin, tup);
                tup[v] = std::min(tup[to], tup[v]);
                auto nbrs = graph.get_all_neighbours(v);
                if (tup[to] > tin[v] &&
                    std::count(nbrs.begin(), nbrs.end(), to) == 1) {
                    bridges.insert({v, to});
                }
            }
        }
    }

    void dfs_find_cutpoints(const Graph::Vertex &v, const Graph::Vertex &ancestor, const Graph &graph,
                            std::vector<vertexColours> &colours, std::set<Graph::Vertex> &cutpoints, size_t &timer,
                            std::vector<int> &tin, std::vector<int> &tup) {
        colours[v] = GRAY;
        tup[v] = tin[v] = timer++;
        int children = 0;
        for (const Graph::Vertex &to : graph.get_all_neighbours(v)) {
            if (to == ancestor) {
                continue;
            } else if (colours[to] == GRAY) {
                tup[v] = std::min(tup[v], tin[to]);
            } else if (colours[to] == WHITE) {
                children++;
                dfs_find_cutpoints(to, v, graph, colours, cutpoints, timer, tin, tup);
                tup[v] = std::min(tup[to], tup[v]);
                if (tup[to] >= tin[v] && ancestor != graph.get_vertex_count()) { // is not root
                    cutpoints.insert(v);
                }
            }
        }
        if (ancestor == graph.get_vertex_count() && children > 1) { // is root
            cutpoints.insert(v);
        }
    }

    void first_dfs_find_strongly_connected_components(const Graph::Vertex &v, const Graph &graph,
                                                      std::vector<vertexColours> &colours,
                                                      std::vector<Graph::Vertex> &order) {
        colours[v] = BLACK;
        for (const Graph::Vertex i : graph.get_all_neighbours(v)) {
            if (colours[i] == WHITE) {
                first_dfs_find_strongly_connected_components(i, graph, colours, order);
            }
        }
        order.push_back(v);
    }

    void second_dfs_find_strongly_connected_components(const Graph::Vertex &v, const Graph &transposed_graph,
                                                       std::vector<vertexColours> &colours,
                                                       std::vector<size_t> &components, size_t cur_component) {
        colours[v] = BLACK;
        components[v] = cur_component;
        for (const Graph::Vertex i : transposed_graph.get_all_neighbours(v)) {
            if (colours[i] == WHITE) {
                second_dfs_find_strongly_connected_components(i, transposed_graph, colours, components, cur_component);
            }
        }
    }


    std::vector<Component> find_components(const Graph &graph) {
        std::vector<Component> connected_components;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                connected_components.emplace_back(std::vector<Graph::Vertex>());
                dfs_find_components(i, graph, colours, connected_components.back());
            }
        }

        return connected_components;
    }

    std::vector<size_t> find_strongly_connected_components(const Graph &graph) {
        std::vector<size_t> components(graph.get_vertex_count());
        size_t cur_component_number = 0;
        std::vector<Graph::Vertex> order;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                first_dfs_find_strongly_connected_components(i, graph, colours, order);
            }
        }

        auto transposed_graph = get_transposed_adj_list_graph(graph);
        std::vector<vertexColours> tr_colours(transposed_graph.get_vertex_count(), WHITE);
        std::reverse(order.begin(), order.end());

        for (Graph::Vertex i = 0; i < transposed_graph.get_vertex_count(); i++) {
            if (tr_colours[order[i]] == WHITE) {
                second_dfs_find_strongly_connected_components(order[i], transposed_graph, tr_colours, components,
                                                              cur_component_number);
                cur_component_number++;
            }
        }
        return components;
    }

    bool is_cicled(const Graph &graph) {
        bool result = 0;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE && !result) {
                result = dfs_check_cicle(i, graph, colours);
            }
        }
        return result;
    }

    std::vector<Graph::Vertex>
    find_shortest_way_bfs(const Graph &graph, const Graph::Vertex &start, const Graph::Vertex &finish) {
        const int NOT_SET = -1;
        std::queue<Graph::Vertex> to_visit;
        std::vector<Graph::Vertex> ancestor(graph.get_vertex_count());
        std::vector<int> distance(graph.get_vertex_count(), NOT_SET);
        to_visit.push(start);
        distance[start] = 0;

        while (!to_visit.empty()) {
            auto cur = to_visit.front();

            to_visit.pop();
            for (const Graph::Vertex &n : graph.get_all_neighbours(cur)) {
                if (distance[n] == NOT_SET) {
                    distance[n] = distance[cur] + 1;
                    ancestor[n] = cur;
                    to_visit.push(n);
                }
            }
        }
        std::vector<Graph::Vertex> result_path;
        auto cur = finish;
        for (int i = 0; i < distance[finish] + 1; i++) {
            result_path.push_back(cur);
            cur = ancestor[cur];
        }
        std::reverse(result_path.begin(), result_path.end());
        return result_path;
    }

    void bfs_find_dist(const Graph &graph, const Graph::Vertex &start, std::vector<size_t> &distances_to_change) {
        const int NOT_SET = -1;
        std::queue<Graph::Vertex> to_visit;
        std::vector<Graph::Vertex> ancestor(graph.get_vertex_count());
        std::vector<int> distance(graph.get_vertex_count(), NOT_SET);
        to_visit.push(start);
        distance[start] = 0;

        while (!to_visit.empty()) {
            auto cur = to_visit.front();

            to_visit.pop();
            for (const Graph::Vertex &n : graph.get_all_neighbours(cur)) {
                if (distance[n] == NOT_SET) {
                    distance[n] = distance[cur] + 1;
                    to_visit.push(n);
                }
            }
        }

        for (int i = 0; i < distance.size(); i++) {
            if (distances_to_change[i] > distance[i]) {
                distances_to_change[i] = distance[i];
            }
        }
    }

    std::vector<Graph::Vertex> graph_top_sort(const Graph &graph) {
        std::vector<Graph::Vertex> top_sort_result;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                dfs_graph_top_sort(i, graph, colours, top_sort_result);
            }
        }

        return top_sort_result;
    }

    std::vector<Graph::Vertex> find_cicle(const Graph &graph) {
        const Graph::Vertex NOT_FOUND = graph.get_vertex_count();
        const Graph::Vertex NOT_SET = graph.get_vertex_count();
        Graph::Vertex last_vertex = NOT_FOUND;
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        std::vector<Graph::Vertex> ancestors(graph.get_vertex_count(), NOT_SET);

        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                last_vertex = dfs_find_cicle(i, graph, colours, ancestors);
            }
            if (last_vertex != NOT_FOUND) {
                break;
            }
        }

        std::vector<Graph::Vertex> cicle;
        if (last_vertex == NOT_FOUND) {
            return cicle;
        }
        cicle.push_back(last_vertex);
        auto cur = ancestors[last_vertex];
        while (cur != last_vertex) {
            cicle.push_back(cur);
            cur = ancestors[cur];
        }
        std::reverse(cicle.begin(), cicle.end());
        return cicle;
    }

    std::set<Edge> find_bridges(const Graph &graph) {
        const int NOT_SET = -1;
        const Graph::Vertex NOT_AC = graph.get_vertex_count();
        std::set<Edge> bridges;
        size_t timer = 0;
        std::vector<int> tin(graph.get_vertex_count(), NOT_SET), tup(graph.get_vertex_count(), NOT_SET);
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                dfs_find_briges(i, NOT_AC, graph, colours, bridges, timer, tin, tup);
            }
        }
        return bridges;
    }

    std::set<Graph::Vertex> find_cutpoints(const Graph &graph) {
        const int NOT_SET = -1;
        const Graph::Vertex NOT_AC = graph.get_vertex_count();
        std::set<Graph::Vertex> cutpoints;
        size_t timer = 0;
        std::vector<int> tin(graph.get_vertex_count(), NOT_SET), tup(graph.get_vertex_count(), NOT_SET);
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);
        for (Graph::Vertex i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                dfs_find_cutpoints(i, NOT_AC, graph, colours, cutpoints, timer, tin, tup);
            }
        }
        return cutpoints;
    }

    bool is_bipartite(const Graph &graph) {
        std::vector<int> part(graph.get_vertex_count(), 0);
        std::vector<vertexColours> colours(graph.get_vertex_count(), WHITE);

        for (int i = 0; i < graph.get_vertex_count(); i++) {
            if (colours[i] == WHITE) {
                part[i] = 1;
                if (!dfs_bipartite(i, graph, colours, part)) {
                    return false;
                }
            }
        }
        return true;
    }

    size_t count_directed_cicle(const Graph &graph) {
        size_t cicle_count = 0;
        // todo
    }

    std::vector<int> ford_find_dist(const WeightedAdjMatrixGraph &graph, const Graph::Vertex &v) {
        const int max_path = 30000;
        std::vector<int> dist(graph.get_vertex_count(), max_path);
        dist[v] = 0;

        for (int k = 0; k < graph.get_edges_count(); ++k) {
            for (int i = 0; i < graph.get_vertex_count(); ++i) {
                for (int u = 0; u < graph.get_vertex_count(); ++u) {
                    if(dist[i] != max_path) {

                        dist[u] = std::min(dist[u], dist[i] + graph.get_weight(i, u));
                    }
                }
            }
        }
        return dist;
    }
};


size_t encode(size_t x, size_t y, size_t N) {
    return x + y * N;
}

std::pair<size_t, size_t> decode(size_t coord, size_t N) {
    return {coord % N, coord / N};
};


int main() {
    size_t vertex_count, edges_count;
    std::cin >> vertex_count >> edges_count;

    WeightedAdjMatrixGraph graph(vertex_count, true);

    for (int i = 0; i < edges_count; i++) {
        Graph::Vertex start, end;
        int w;
        std::cin >> start >> end >> w;
        graph.add_weighted_edge(start - 1, end - 1, w);
    }

    for (int i : GraphAlgorithms::ford_find_dist(graph, 0)) {
        std::cout << i << ' ';
    }
}