/*******************************************
 * Author: Shuai Zhou.
 * All Rights Reserved.
 *******************************************/
#include "mapfaa_ppsipp.hpp"


namespace raplab{

    // Prioritize Planning using SIPP (Safe Interval Path Planning) as the single agent planner
    PP_SIPP::PP_SIPP() {}

    PP_SIPP::~PP_SIPP() {
        _graph = nullptr;
        _starts.clear();
        _goals.clear();
        _paths.clear();
    }

    CostVec PP_SIPP::GetPlanCost(long nid) {
        CostVec out(_graph->CostDim(), 0);
        return out;
    }

    TimePathSet PP_SIPP::GetPlan(long nid) {
        return _paths;
    }

    std::unordered_map<std::string, double> PP_SIPP::GetStats() {
        std::unordered_map<std::string, double> stats;
        return stats;
    }

    int PP_SIPP::Solve(std::vector<long>& starts, std::vector<long>& goals, double time_limit, double eps) {
        _starts = starts;
        _goals = goals;
        _time_limit = time_limit;
        _cost = 0;
        return 1;
    }

    std::vector<std::unordered_map<long, double>> PP_SIPP::generate_distable() {
        std::vector<std::unordered_map<long, double>> distable;
        for (size_t i = 0; i < _starts.size(); ++i) {
            distable.push_back(generate_single_dis_table(static_cast<int>(i)));
        }
        return distable;
    }

    std::unordered_map<long, double> PP_SIPP::generate_single_dis_table(int agent) {
        std::deque<long> tmp;
        tmp.push_back(_starts[agent]);
        const double inf = std::numeric_limits<double>::max();
        std::unordered_map<long, double> dist_table;

        dist_table[_goals[agent]] = 0;

        while (!tmp.empty()) {
            long curr = tmp.front();
            tmp.pop_front();
            std::vector<long> successors = _graph->GetSuccs(curr);
            for (long neigh : successors) {
                if (!edge_cost.empty() && edge_cost.find(agent) != edge_cost.end() && edge_cost[agent].find(edge_hash(curr, neigh)) != edge_cost[agent].end()) {
                    // this edge cost for this agent is specified
                    if (dist_table.find(neigh) == dist_table.end() || dist_table[neigh] > dist_table[curr] + edge_cost[agent].at(edge_hash(curr, neigh))) {
                        dist_table[neigh] = dist_table[curr] + edge_cost[agent].at(edge_hash(curr, neigh));
                        tmp.push_back(neigh);
                    }
                } else {
                    // this edge cost for this agent is not specified so using the default duration cost
                    if (dist_table.find(neigh) == dist_table.end() || dist_table[neigh] > dist_table[curr] + _duration[agent]) {
                        dist_table[neigh] = dist_table[curr] + _duration[agent];
                        tmp.push_back(neigh);
                    }
                }
            }
        }
        return dist_table;
    }

    void PP_SIPP::set_duration(std::vector<double> duration) {
        _duration = duration; //duration of each agent
    }

    double PP_SIPP::get_h(int agent, long vertex) {
        return _dis_table[agent].at(vertex)* _duration[agent];
    }

    double PP_SIPP::get_D(int agent, long v_departure, long v_arrive) {
        auto edge = edge_hash(v_departure,v_arrive);
        if (!edge_cost.empty() && edge_cost.find(agent) != edge_cost.end() && edge_cost.at(agent).find(edge) != edge_cost.at(agent).end())
        {
            return edge_cost.at(agent).at(edge);
        }
        return _duration[agent];
    }

    bool PP_SIPP::run_PP() {
        // basic prioritize planning
        _dis_table = generate_distable();
        _start_time = std::chrono::steady_clock::now();
        for (int i = 0; i < _starts.size(); ++i) {
            auto success = run_SIPP(i);
            if (!success) {
                return false; // If any agent fails, return false
            }
        }
        _runtime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - _start_time).count();
        // If all agents successfully find a path, return true
        std::cout<<"Path Found! Sum Of Cost: "<<_cost<<" Makespan: "<<_makespan<<std::endl;
        return true;
    }

    bool PP_SIPP::run_SIPP(int agent) {
        auto safe_intervals = get_safe_intervals(_starts[agent]);
        if (safe_intervals.empty() || safe_intervals[0].first > 0) {
            return false; // No safe intervals available or the first interval starts after time 0, that means the agent cannot start
        }
        auto n_0 = new Sipp_Node(_starts[agent],0.0, safe_intervals[0], get_h(agent, _starts[agent]));
        Open.clear();
        Closed.clear();
        Visited.clear();
        // closed set distinguish between the nodes by their vertex _v and arrival time _t
        Open.insert(n_0);

        while (!Open.empty()) {
            Sipp_Node* n = *Open.begin();
            Open.erase(Open.begin());
            Closed.insert(*n);

            // Check if the current node is the goal
            if (n->_v == _goals[agent] && n->_interval.second == std::numeric_limits<double>::infinity()) {
                // which means it can forever stay at goal
                // Reconstruct the path and return true
                reconstruct_Path(agent, n);
                _cost += n->_t;
                if (n->_t > _makespan) {
                    _makespan = n->_t;
                }
                return true;
            }

            if (Timeout()) {
                std::cout << "Timeout!" << std::endl;
                return false; // Timeout
            }

            // Expand the node and generate successors

            auto successors = get_successors(n, agent);
            if (successors.empty()) {
                continue; // No successors available
            }
            for (const auto& n_new : successors) {
                auto pair = std::make_pair(n_new._v, n_new._interval);
                if (Visited.find(pair) == Visited.end()) {
                    Visited[pair] = n_new;
                    Open.insert(new Sipp_Node(n_new)); // Push the new node to the open set
                } else {
                    // Check if the new node is better than the existing one
                    auto existing_node = Visited[pair];
                    if (n_new._t < existing_node._t) {
                        // Update the existing node with the new one
                        remove_node_from_open(&existing_node);
                        Visited[pair] = n_new;
                        Open.insert(new Sipp_Node(n_new)); // Push the new node to the open set
                    }
                }
            }
        }

        return false;
    }

    std::vector<std::pair<double, double>> PP_SIPP::get_safe_intervals(long vertex) {
        if (_unsafe_intervals.find(vertex) != _unsafe_intervals.end()) {
            std::vector<std::pair<double, double>> safe_intervals;
            // Find the gaps between unsafe intervals
            double last_end = 0.0;
            for (const auto& interval : _unsafe_intervals[vertex]) {
                if (interval.first > last_end) {
                    safe_intervals.emplace_back(last_end, interval.first);
                }
                last_end = std::max(last_end, interval.second);
            }
            // Add the last safe interval
            if (last_end != std::numeric_limits<double>::infinity()) {
                safe_intervals.emplace_back(last_end, std::numeric_limits<double>::infinity());
            }
            return safe_intervals;
        } else {
            // If there are no unsafe intervals for this vertex, return a single safe interval
            return {{0.0, std::numeric_limits<double>::infinity()}};
        }
    }

    std::vector<Sipp_Node> PP_SIPP::get_successors(Sipp_Node *n, int agent) {
        std::vector<Sipp_Node> successors;
        for (const auto& m: _graph->GetSuccs(n->_v)) {
            double m_time = get_D(agent,n->_v,m);
            // for mapfaa, we consider whether the time between the end of safe interval and _t
            // is enough for the agent to move to the next node
            if (n->_interval.second - n->_t < m_time) {
                continue; // Not enough time to move to the next node
            }
            // for mapfaa the end_t, it should be the end of the safe interval
            // it means the action should be start after _t and end before the end of the safe interval
            double start_t = n->_t + m_time;
            double end_t = n->_interval.second;
            auto safe_intervals = get_safe_intervals(m);
            if (safe_intervals.empty()) {continue;} // No safe intervals available
            // Check if the new interval overlaps with any unsafe intervals
            for (const auto& interval : safe_intervals) {
                if (interval.first > end_t || interval.second < start_t) {
                    continue; // No overlap with the safe interval
                }
                // just to generate a search node first  where t is not the real t we want to set
                Sipp_Node n_new(m, 0.0, interval, 0.0, n);
                if (Closed.find(n_new) != Closed.end()) {
                    continue; // Already in closed set
                }
                // for mapf_aa, if safe interval of m smaller than m_time,
                // we can not move to the next node
                if (interval.second - interval.first < m_time) {continue;} // Not enough time to move to the next node
                // if t exists
                // which means the beginning of safe interval smaller than latest departure time
                // and the end of safe interval larger than the earliest arrival time
                if (interval.first <= end_t - m_time && interval.second >= start_t) {
                    // the earliest arrival time
                    double t = std::max(n->_t, interval.first) + m_time;
                    // Create a new Sipp_Node for the successor
                    n_new._t = t;
                    n_new._f = t + get_h(agent, m);
                    // Add the new node to the successors
                    successors.push_back(n_new);
                } else {
                    // no such t exists
                    continue;
                }
            }
        }
        return successors;
    }


    void PP_SIPP::reconstruct_Path(int agent, raplab::Sipp_Node *n) {
        std::vector<std::tuple<long, long, double, double>> single_path;
        // last n.v should last to infinity
        auto goal_wait = std::make_pair(n->_t, std::numeric_limits<double>::infinity());
        add_to_unsafe_intervals(n->_v, goal_wait);
        while (n->_parent != nullptr) {
            double m_time = get_D(agent,n->_v,n->_parent->_v);
            double wait_time = n->_t - n->_parent->_t - m_time;
            if (wait_time == 0) {
                single_path.emplace_back(n->_parent->_v, n->_v, n->_parent->_t, n->_parent->_t + m_time);
                // add the unsafe interval to the unsafe intervals
                auto move = std::make_pair(n->_parent->_t, n->_parent->_t + m_time);
                add_to_unsafe_intervals(n->_parent->_v, move);
                add_to_unsafe_intervals(n->_v, move);
            } else {
                single_path.emplace_back(n->_parent->_v, n->_v, n->_t - m_time, n->_t);
                single_path.emplace_back(n->_parent->_v, n->_parent->_v,n->_parent->_t,n->_parent->_t + wait_time);
                // add the unsafe interval to the unsafe intervals
                auto move = std::make_pair(n->_t - m_time, n->_t);
                add_to_unsafe_intervals(n->_v, move);
                auto waitandmove = std::make_pair(n->_parent->_t, n->_t);
                add_to_unsafe_intervals(n->_parent->_v, waitandmove);
            }
            n = n->_parent;
        }
        // reverse the path
        std::reverse(single_path.begin(), single_path.end());
        // add the path to the final path
        _all_paths.emplace_back(single_path);
    }
}
