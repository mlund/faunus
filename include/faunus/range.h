#ifndef CONTINUOUSRANGE_H
#define CONTINUOUSRANGE_H

/** 
 * @brief Iterable class for a continuous range of integers
 *
 * Example:
 *
 *     ContinuousRange<int> r(2,5); // first, last 
 *     for (auto i : r)             // iterator access
 *       std::cout << i;            // -> 2345
 *     r.front();                   // -> 2
 *     r.back();                    // -> 5
 *     r.size();                    // -> 4
 *     r.resize( r.size()+1 );      // -> size=5, back=6.
 *
 * For more information see:
 *
 * - <http://stackoverflow.com/questions/7185437/is-there-a-range-class-in-c0x-aka-c11-for-use-with-range-based-for-loops>
 * - <https://bitbucket.org/AraK/range/wiki/Home>
 */
template<class T=int>
class ContinuousRange {
  public:
    class iterator : public std::iterator<std::forward_iterator_tag,T> {
      private:
        friend class ContinuousRange;
        T i_;
      public:

        T operator*() const { return i_; }

        const iterator &operator++()
        {
          ++i_;
          return *this;
        }

        iterator operator++(int)
        {
          iterator copy=*this;
          ++i_;
          return copy;
        }

        iterator operator+(T i)
        {
          iterator copy=*this;
          copy.i_+=i;
          return copy;
        }

        iterator operator-(T i)
        {
          iterator copy=*this;
          copy.i_-=i;
          return copy;
        }

        bool operator==(const iterator &other) const { return i_ == other.i_; }
        bool operator!=(const iterator &other) const { return i_ != other.i_; }

      protected:
        iterator(T start=0) : i_(start) {}
    };

    iterator begin() const { return begin_; }    //!< Iterator to beginning

    iterator end() const { return end_; }        //!< Iterator to end (end points to last element + 1)

    T front() const { return begin_.i_; }        //!< Get first value in range

    T back() const { return end_.i_-1; }         //!< Get last value in range

    bool empty() const { return (end_.i_<=begin_.i_) ? true : false; } //!< Determines if range is empty.

    void clear() { *this = ContinuousRange<T>(); } //!< Clear range

    /** @brief Resize range, keeping same beginning */
    void resize(T newsize)
    {
      assert( newsize >= 0 );
      if ( front()<0 )
        begin_.i_ = 0;
      end_.i_ = begin_.i_ + newsize;

      assert( size() == newsize );
    }

    T size() const { return end_.i_-begin_.i_; } //!< Size of range

    void setfront(T front) { begin_.i_=front; }  //!< Set first element

    void setback(T back) { end_.i_=back+1; }     //!< Set last element

    /**
     * @brief Set range `[front:back]`
     * @param front First particle
     * @param back Last particle
     */
    void setrange(T front, T back=-1) 
    {
      setfront(front);
      if (back>=0)
        setback(back);
      else setback(front-1);
      assert(size()>=0);
    }

    bool find(T i) const                         //!< Check if index is part of group
    {
      if (i>back()) return false;
      if (i<front()) return false;
      return true;
    }

    ContinuousRange(T first=0, T size=0) : begin_(first), end_(first+size) {
      static_assert( std::is_integral<T>::value, "T must be of integral type" );
    }

  private:

    iterator begin_;
    iterator end_;
};
#endif
