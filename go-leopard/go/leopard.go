package main

import "C"

import (
	"unsafe"

	"github.com/klauspost/reedsolomon"
)

//export encode
func encode(ptr **byte, shard_size uint, data_shards uint, total_shards uint) int {
	input := unsafe.Slice(ptr, total_shards)
	shards := make([][]byte, total_shards)

	// copy shards to go
	for i := uint(0); i < total_shards; i++ {
		shard := unsafe.Slice(input[i], shard_size)
		shards[i] = make([]byte, shard_size)
		copy(shards[i], shard)
	}

	// create the codec and encode
	codec, err := reedsolomon.New(int(data_shards), int(total_shards-data_shards), reedsolomon.WithLeopardGF(true))
	if err != nil {
		return -1
	}
	if codec.Encode(shards) != nil {
		return -2
	}

	// copy shards back to the input ptr
	for i := uint(0); i < total_shards; i++ {
		shard := unsafe.Slice(input[i], shard_size)
		copy(shard, shards[i])
	}

	return 0
}

func main() {}
