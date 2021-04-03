#import "Basic";




main :: () {
    // data : []float64;
    data := NewArray(4, float64, initialized=false);
    defer free(data.data);
    print("%\n", data);
}