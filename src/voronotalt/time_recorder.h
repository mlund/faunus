#ifndef VORONOTALT_TIME_RECORDER_H_
#define VORONOTALT_TIME_RECORDER_H_

namespace voronotalt
{

class TimeRecorder
{
public:
	TimeRecorder()
	{
	}

	virtual ~TimeRecorder() noexcept
	{
	}

	virtual void reset() noexcept
	{
	}

	virtual void record_elapsed_miliseconds_and_reset(const char* /*message*/) noexcept
	{
	}
};

}

#endif /* VORONOTALT_TIME_RECORDER_H_ */
