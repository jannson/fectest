package main

import (
	"os"

	"fmt"
	"github.com/klauspost/reedsolomon"
)

func main() {
	dataShards := 11
	parShards := 6
	enc, err := reedsolomon.New(dataShards, parShards)
	checkErr(err)

	data := make([]byte, dataShards)
	copy(data, []byte("hello world"))

	shards, err := enc.Split(data)
	checkErr(err)
	fmt.Println(shards)

	err = enc.Encode(shards)
	checkErr(err)
	fmt.Println(shards)

	dec, err := reedsolomon.New(dataShards, parShards)
	checkErr(err)
	//fmt.Println(len(shards), len(shards[0]))

	shards[1] = nil
	shards[3] = nil
	shards[4] = nil
	shards[5] = nil
	//fmt.Println(shards)

	err = dec.Reconstruct(shards)
	checkErr(err)
	fmt.Println(shards)
}

func checkErr(err error) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err.Error())
		os.Exit(2)
	}
}
