#include "DeltaSteppingSequential.h"

DeltaSteppingSequential::DeltaSteppingSequential(const Graph& graph, const int source, const bool is_verbose)
    : graph_(graph), source_(source), is_verbose_(is_verbose) {
    
    delta_ = 1.0 / (static_cast<double>(graph_.getMaxDegree()) + 1);
    dist_.resize(graph_.getNumVertices(), std::numeric_limits<double>::infinity());
    pred_.resize(graph_.getNumVertices(), -1);

    light_edges_.resize(graph_.getNumVertices());
    heavy_edges_.resize(graph_.getNumVertices());

    computeLightAndHeavyEdges();

    int bucket_size = static_cast<int>(graph_.getNumVertices() / delta_) + 1;
    buckets_.resize(bucket_size);

    if (is_verbose_) {
        std::cout << "Bucket size: " << bucket_size << std::endl;
    }

    for (int i = 0; i < graph_.getNumVertices(); i++) {
        if (i != source_) {
            buckets_[bucket_size - 1].insert(i);
        }
    }

    buckets_[0].insert(source_);
    dist_[source_] = 0;
    pred_[source_] = source_;

    if (is_verbose_) {
        printLightAndHeavyEdges();
        printAllBuckets();
    }
}

void DeltaSteppingSequential::computeLightAndHeavyEdges() {
    for (const Edge& edge : graph_.getEdges()) {
        if (edge.getWeight() <= delta_) {
            light_edges_[edge.getFrom()].push_back(edge.getTo());
        } else {
            heavy_edges_[edge.getFrom()].push_back(edge.getTo());
        }
    }
}

void DeltaSteppingSequential::findBucketRequests(const std::set<int>& bucket, std::vector<Edge>* light_requests, std::vector<Edge>* heavy_requests) {
    for (int vertex_id : bucket) {
        buckets_[bucket_counter_].erase(vertex_id);
        if (is_verbose_) {
            std::cout << "Erased " << vertex_id << " from bucket " << bucket_counter_ << std::endl;
            printBucket(bucket_counter_);
        }

        for (int l_edge_vertex_id : light_edges_[vertex_id]) {
            light_requests->emplace_back(vertex_id, l_edge_vertex_id, graph_.getEdgeWeight(vertex_id, l_edge_vertex_id));
        }

        for (int h_edge_vertex_id : heavy_edges_[vertex_id]) {
            heavy_requests->emplace_back(vertex_id, h_edge_vertex_id, graph_.getEdgeWeight(vertex_id, h_edge_vertex_id));
        }
    }
}

void DeltaSteppingSequential::relax(const Edge& selected_edge) {
    int from_vertex = selected_edge.getFrom();
    int to_vertex = selected_edge.getTo();
    double edge_weight = selected_edge.getWeight();
    double tentative_dist = dist_[from_vertex] + edge_weight;

    if (tentative_dist < dist_[to_vertex]) {
        int i = static_cast<int>(std::floor(dist_[to_vertex] / delta_));
        int j = static_cast<int>(std::floor(tentative_dist / delta_));

        if (i < buckets_.size() && i >= 0) {
            buckets_[i].erase(to_vertex);
        }
        if (j < buckets_.size() && j >= 0) {
            buckets_[j].insert(to_vertex);
        }

        dist_[to_vertex] = tentative_dist;
        pred_[to_vertex] = from_vertex;
    }
}

void DeltaSteppingSequential::resolveRequests(std::vector<Edge>* requests) {
    for (const Edge& request : *requests) {
        relax(request);
    }
}

void DeltaSteppingSequential::solve() {
    while (bucket_counter_ < buckets_.size()) {
        std::set<int> current_bucket = buckets_[bucket_counter_];
        while (!current_bucket.empty()) {
            std::set<int> current_bucket_update;

            for (int vertex_id : current_bucket) {
                for (int neighbor_vertex : light_edges_[vertex_id]) {
                    relax(Edge(vertex_id, neighbor_vertex, graph_.getEdgeWeight(vertex_id, neighbor_vertex)));
                    current_bucket_update.insert(neighbor_vertex);
                }
            }

            current_bucket = current_bucket_update;
        }
        bucket_counter_++;
        while (bucket_counter_ < buckets_.size() && buckets_[bucket_counter_].empty()) {
            bucket_counter_++;
        }
    }
}

void DeltaSteppingSequential::solveLightHeavy() {
    while (bucket_counter_ < buckets_.size()) {
        std::vector<Edge> light_requests, heavy_requests;
        std::set<int> current_bucket = buckets_[bucket_counter_];
        while (!current_bucket.empty()) {
            findBucketRequests(current_bucket, &light_requests, &heavy_requests);

            resolveRequests(&light_requests);
            light_requests.clear();

            current_bucket = buckets_[bucket_counter_];
        }

        resolveRequests(&heavy_requests);
        heavy_requests.clear();

        bucket_counter_++;
        while (bucket_counter_ < buckets_.size() && buckets_[bucket_counter_].empty()) {
            bucket_counter_++;
        }
    }
}

void DeltaSteppingSequential::printSolution() const {
    std::cout << "Solution: " << std::endl;
    for (int i = 0; i < graph_.getNumVertices(); i++) {
        std::cout << "Vertex " << i << ": " << dist_[i] << std::endl;
    }
}

void DeltaSteppingSequential::printLightAndHeavyEdges() const {
    std::cout << "Light edges: " << std::endl;
    for (int i = 0; i < graph_.getNumVertices(); i++) {
        std::cout << "Vertex " << i << ": ";
        for (int light_edge : light_edges_[i]) {
            std::cout << light_edge << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Heavy edges: " << std::endl;
    for (int i = 0; i < graph_.getNumVertices(); i++) {
        std::cout << "Vertex " << i << ": ";
        for (int heavy_edge : heavy_edges_[i]) {
            std::cout << heavy_edge << " ";
        }
        std::cout << std::endl;
    }
}

void DeltaSteppingSequential::printAllBuckets() const {
    for (size_t bucket_id = 0; bucket_id < buckets_.size(); bucket_id++) {
        printBucket(bucket_id);
    }
}

void DeltaSteppingSequential::printBucket(size_t bucket_id) const {
    std::cout << "Bucket [" << bucket_id << "], size " << buckets_[bucket_id].size() << ": ";
    for (int bucket_item : buckets_[bucket_id]) {
        std::cout << bucket_item << " ";
    }
    std::cout << std::endl;
}
