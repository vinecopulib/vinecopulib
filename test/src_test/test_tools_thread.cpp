// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <atomic>
#include <chrono>
#include <stdexcept>
#include <thread>
#include <vector>
#include <vinecopulib/misc/tools_interface.hpp>

namespace test_tools_thread {

using namespace vinecopulib::tools_thread;

// The pool's state is only consistent to the thread holding its mutex, so the
// jobs below touch nothing but an atomic counter and every assertion is made
// from the main thread after a wait(), join() or the destructor.

TEST(test_tools_thread, runs_every_job)
{
  const int n_jobs = 1000;
  const std::vector<size_t> worker_counts{ 0, 1, 4 };
  for (size_t n_workers : worker_counts) {
    std::atomic<int> ran{ 0 };
    ThreadPool pool(n_workers);
    for (int i = 0; i < n_jobs; ++i)
      pool.push([&ran]() noexcept { ++ran; });
    pool.wait();
    EXPECT_EQ(ran.load(), n_jobs) << "with " << n_workers << " workers";
    pool.join();
  }
}

// Without workers the job runs at the push site, and so does its exception.
TEST(test_tools_thread, no_workers_runs_inline)
{
  ThreadPool pool(0);
  std::thread::id job_id;
  pool.push([&job_id]() noexcept { job_id = std::this_thread::get_id(); });
  EXPECT_EQ(job_id, std::this_thread::get_id());

  EXPECT_THROW(pool.push([] { throw std::runtime_error("job failed"); }),
               std::runtime_error);
  EXPECT_NO_THROW(pool.wait());
  EXPECT_NO_THROW(pool.join());
}

TEST(test_tools_thread, exception_propagates_from_wait)
{
  ThreadPool pool(2);
  pool.push([] { throw std::runtime_error("job failed"); });
  try {
    pool.wait();
    FAIL() << "wait() did not rethrow the job's exception";
  } catch (const std::runtime_error& e) {
    EXPECT_STREQ(e.what(), "job failed");
  }
}

// An error cancels what has not started, but wait() still returns only once the
// jobs already running have finished.
TEST(test_tools_thread, error_cancels_queued_jobs)
{
  const int n_jobs = 200;
  std::atomic<int> ran{ 0 };
  int ran_on_return = 0;
  {
    ThreadPool pool(2);
    pool.push([] { throw std::runtime_error("job failed"); });
    for (int i = 0; i < n_jobs; ++i)
      pool.push([&ran] {
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        ++ran;
      });
    EXPECT_THROW(pool.wait(), std::runtime_error);
    ran_on_return = ran.load();
    EXPECT_LT(ran_on_return, n_jobs);
  }
  EXPECT_EQ(ran.load(), ran_on_return) << "a job ran after wait() returned";
}

// A worker may only exit once the queue is drained, so the destructor runs
// every job that was pushed.
TEST(test_tools_thread, destructor_drains_queued_jobs)
{
  const int n_jobs = 100;
  std::atomic<int> ran{ 0 };
  {
    ThreadPool pool(2);
    for (int i = 0; i < n_jobs; ++i)
      pool.push([&ran]() noexcept { ++ran; });
  }
  EXPECT_EQ(ran.load(), n_jobs);
}

TEST(test_tools_thread, clear_then_join_terminates)
{
  std::atomic<int> ran{ 0 };
  ThreadPool pool(2);
  for (int i = 0; i < 1000; ++i)
    pool.push([&ran]() noexcept { ++ran; });
  pool.clear();
  EXPECT_NO_THROW(pool.join());
  EXPECT_LE(ran.load(), 1000);
}

TEST(test_tools_thread, push_after_join_throws)
{
  ThreadPool pool(2);
  pool.join();
  EXPECT_THROW(pool.push([]() noexcept {}), std::runtime_error);
  EXPECT_NO_THROW(pool.join());
}

// How VinecopSelector drives its pool: map and wait, once per tree.
TEST(test_tools_thread, reusable_across_rounds)
{
  const int n_items = 100;
  std::atomic<int> ran{ 0 };
  ThreadPool pool(4);
  for (int tree = 0; tree < 5; ++tree) {
    std::vector<int> items(n_items, 0);
    pool.map([&ran](int) noexcept { ++ran; }, items);
    pool.wait();
    EXPECT_EQ(ran.load(), n_items * (tree + 1));
  }
  pool.join();
}

// Repeated start/stop cycles, to be run under -fsanitize=thread.
TEST(test_tools_thread, repeated_shutdown_is_consistent)
{
  const int n_jobs = 50;
  for (int rep = 0; rep < 50; ++rep) {
    std::atomic<int> ran{ 0 };
    ThreadPool pool(4);
    for (int i = 0; i < n_jobs; ++i)
      pool.push([&ran]() noexcept { ++ran; });
    pool.join();
    EXPECT_EQ(ran.load(), n_jobs);
  }
}
}
