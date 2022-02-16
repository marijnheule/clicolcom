
#include <iostream>
#include "interval_list.hpp"



gc::interval_list::interval_list() { 
		inf.push_back(-0x1000000);
		inf.push_back(0x1000000);
		sup.push_back(-0x1000000);
		sup.push_back(0x1000000);
		next.push_back(1);
		next.push_back(-1);
	
		size_ = 0; 
}

bool gc::interval_list::contain(const int x)
{
    int pos = next[0];
    while (pos > 1 && sup[pos] < x) {
        pos = next[pos];
    }
    return inf[pos] <= x;
}

bool gc::interval_list::remove(const int x) {
	int prev=0, pos=next[0], idx;
	while(pos > 1 && sup[pos] < x) {
			prev = pos;
			pos = next[pos];
	}
	if(inf[pos] > x) {	
		  if(sup[prev] == inf[pos] - 2) { // merge
				sup[prev] = sup[pos];
				next[prev] = next[pos];
				freed.push_back(pos);
			} else if(sup[prev] == x - 1) {
				sup[prev] = x;
			} else if(inf[pos] == x + 1) {
				inf[pos] = x;
			} else { // new singleton interval
				if(freed.size() == 0) {
						inf.push_back(x);
						sup.push_back(x);
						idx = next.size();
						next.push_back(pos);
				} else {
						idx = freed.back();
						freed.pop_back();
						inf[idx] = x;
						sup[idx] = x; 
						next[idx] = pos;
						next[prev] = idx;
				}
				next[prev] = idx;
			}
		
			++size_;
			return true;
	}
	
	return false;
}

// bool gc::interval_list::remove(const int x) {
// 	int prev=0, pos=next[0], idx;
// 	while(pos > 1 && sup[pos] < x) {
// 			prev = pos;
// 			pos = next[pos];
// 	}
// 	if(inf[pos] > x) {
// 		  if(sup[prev] == inf[pos] - 2) { // merge
// 				sup[prev] = sup[pos];
// 				next[prev] = next[pos];
// 				freed.push_back(pos);
// 			} else if(sup[prev] == x - 1) {
// 				sup[prev] = x;
// 			} else if(inf[pos] == x + 1) {
// 				inf[pos] = x;
// 			} else { // new singleton interval
// 				if(freed.size() == 0) {
// 						inf.push_back(x);
// 						sup.push_back(x);
// 						idx = next.size();
// 						next.push_back(pos);
// 				} else {
// 						idx = freed.back();
// 						freed.pop_back();
// 						inf[idx] = x;
// 						sup[idx] = x;
// 						next[idx] = pos;
// 						next[prev] = idx;
// 				}
// 				next[prev] = idx;
// 			}
//
// 			++size;
// 			return true;
// 	}
//
// 	return false;
// }


int gc::interval_list::min() const {
		int c = sup[next[0]];
		return (c == 0x1000000 ? 0 : c + 1); 
}

void gc::interval_list::fill() {
	inf.resize(2);
	sup.resize(2);
	next.resize(2);
	next[0] = 1;
	next[1] = -1;
	
	size_ = 0;
  freed.clear();
	
	check_consistency();
}

/*!@name Printing*/
//@{
std::ostream& gc::interval_list::display(std::ostream& os) const
{
	
		os << "{";
		int pos=next[0];
		while(pos > 1) {
				if(inf[pos] == sup[pos])
						os << " " << inf[pos];
				else
					os << " [" << inf[pos] << ".." << sup[pos] << "]";
				pos = next[pos];
		}
	
    os << " }";
		
		// os << std::endl;
		// for( auto s : inf ) {
		// 	os << std::setw(2) << s << " ";
		// }
		// os << std::endl;
		// for( auto s : sup ) {
		// 	os << std::setw(2) << s << " ";
		// }
		// os << std::endl;
		// for( auto s : next ) {
		// 	os << std::setw(2) << s << " ";
		// }
		// os << std::endl;
		
		
    return os;
}
//@}

void gc::interval_list::check_consistency() {
	
	assert( inf[0] == -0x1000000 ) ;
	assert( inf[1] == 0x1000000 ) ;
	assert( sup[0] == -0x1000000 ) ;
	assert( inf[1] == 0x1000000 ) ;
	assert( next[0] == 1 ) ;
	assert( next[1] == -1 ) ;
	
	int elt{0}, nxt{next[0]};
	size_t count{0};
	while(nxt != -1) {
		assert(inf[elt] <= sup[elt]);
		// int nxt{next[elt]};
		assert(sup[elt]+1 < inf[nxt]);
		// elt = nxt;
		++count;
		
		elt = nxt;
		nxt = next[nxt];
	}
	
	assert(count == size_ + 1);
	
}


std::ostream& gc::operator<<(std::ostream& os, const gc::interval_list& x)
{
    return x.display(os);
}

std::ostream& gc::operator<<(std::ostream& os, const gc::interval_list* x)
{
    return (x ? x->display(os) : os);
}
