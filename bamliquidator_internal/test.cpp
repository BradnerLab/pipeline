#include <iostream>
#include <iomanip>
#include <deque>
#include <vector>
#include <stdexcept>

struct ReadItem
{
  unsigned int start;
  /* read stop is start + strlen(seq)
  this *stop* will only be used for computing density
  will not be reported to js for bed plotting
  the actual stop need to be determined by cigar
  */
  unsigned int stop;
  uint32_t flag; // flag from bam
  char strand;
  std::vector<uint32_t> cigar;
};

int intMin(int a, int b)
{
  if(a < b) return a;
  return b;
}

int intMax(int a, int b)
{
  if(a > b) return a;
  return b;
}   

std::deque<ReadItem> sample(unsigned int start, unsigned int stop)
{
  std::deque<ReadItem> items;

  for (size_t i=start; i<stop; ++i)
  {
    ReadItem r;
    r.start = i;
    r.stop = i + 1;
    items.push_back(r);
  }

  return items;
}

template <typename T>
void print(T* startArr, T* stopArr, std::vector<double> data, std::deque<ReadItem> items, size_t spnum)
{
  /*std::cout << "Read items: " << std::endl;
  for (ReadItem item : items)
  {
    std::cout << "[" << item.start << ", " << item.stop << "): " << 1 << std::endl;
  }
  */

  double sum = 0;
  std::cout << std::endl << "Density values: " << std::endl;
  for (size_t i=0; i < spnum; ++i)
  {
    std::cout << "[" << startArr[i] << ", " << stopArr[i] << "): " << data[i] << std::endl;
    sum += data[i];
  }
  std::cout << "Sum: " << sum << std::endl;

  if (spnum != data.size())
  {
    throw std::runtime_error("something isn't right");
  }
}

const unsigned int start=76619699, stop=76620366, spnum=200;

void test_current_double_version()
{
  std::cout << std::endl << "current double version:" << std::endl;
  std::vector<double> data(spnum, 0);
  std::deque<ReadItem> items = sample(start, stop);

  double startArr[spnum], stopArr[spnum];
  double pieceLength = (double)(stop-start) / (double)spnum;
  for(int i=0; i<spnum; i++)
  {
    startArr[i] = (double)(start + pieceLength*i);
    stopArr[i] = (double)(start + pieceLength*(i+1));
  }

  for(const ReadItem& item : items)
  {
    // collapse this bed item onto the density counter
    for(int i=0; i<spnum; i++)
    {
      if(item.start > stopArr[i]) continue;
      if(item.stop < startArr[i]) break;
      int start=intMax(item.start,(int)startArr[i]);
      int stop=intMin(item.stop,(int)stopArr[i]);
      if(start<stop)
      {
        // as Charles suggested, add the fraction of the read (overlapping with the bin)
        // instead of just counting the read
        data[i] += stop-start;
      }
    }
  }

  print(startArr, stopArr, data, items, spnum);
}

void test_proposed_int_version()
{
  std::cout << std::endl << "proposed int version:" << std::endl;
  std::vector<double> data(spnum, 0);
  std::deque<ReadItem> items = sample(start, stop);

  int startArr[spnum], stopArr[spnum];
  int pieceLength = (int)(stop-start) / spnum;
  for(int i=0; i<spnum; i++)
  {
    startArr[i] = (int)(start + pieceLength*i);
    stopArr[i] = (int)(start + pieceLength*(i+1));
  }

  for(const ReadItem& item : items)
  {
    // collapse this bed item onto the density counter
    for(int i=0; i<spnum; i++)
    {
      if(item.start > stopArr[i]) continue;
      if(item.stop < startArr[i]) break;
      int start=intMax(item.start,(int)startArr[i]);
      int stop=intMin(item.stop,(int)stopArr[i]);
      if(start<stop)
      {
        // as Charles suggested, add the fraction of the read (overlapping with the bin)
        // instead of just counting the read
        data[i] += stop-start;
      }
    }
  }

  print(startArr, stopArr, data, items, spnum);
}

int main()
{
  std::cout << std::fixed << std::setprecision(2); 
  test_current_double_version();
  //test_proposed_int_version();

  return 0;
}
