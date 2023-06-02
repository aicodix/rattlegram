#include "encoder.hh"
#include <cstdint>

static EncoderInterface *encoder;
//static DecoderInterface *decoder;
typedef DSP::Complex<float> cmplx;
int RATE = 8000;
int SYMBOL_LENGHT = (1280 * RATE) / 8000;
int GUARD_LENGTH = SYMBOL_LENGHT/8;
int EXTENDED_LENGTH = SYMBOL_LENGHT + GUARD_LENGTH;
int CHANNEL = 0;
int DEFAULT_NOISE_SYMBOLS = 1;
const char *CALL_SIGN = "ANONYMOUS";
const char *MESSAGE = "TEST";
int CARRIER_FREQUENCY = 1500;
const char *OUTPUT_FILE = "test.wav";
int OUTPUT_RATE = 8000;
int OUTPUT_BITS = 16;

int main(int argc, char **argv) {
    DSP::WriteWAV<float> output_file(OUTPUT_FILE, OUTPUT_RATE, OUTPUT_BITS, 1);
    encoder = new(std::nothrow) Encoder<8000>(&output_file);
    int16_t buffer[10*EXTENDED_LENGTH] = {0};
    encoder->configure(reinterpret_cast<const uint8_t*>(&MESSAGE[0]), 
        reinterpret_cast<const int8_t*>(&CALL_SIGN[0]), 
        CARRIER_FREQUENCY, 
        DEFAULT_NOISE_SYMBOLS, 
        false);
    for (int i = 0; i < 5; ++i) {
        encoder->produce(&buffer[i*2*EXTENDED_LENGTH], 0);
	}
    
    for (int i = 0; i < 10*1440; i++) {
        std::cout << buffer[i] << " ";
    }

    DSP::WriteWAV<int16_t> out(OUTPUT_FILE, OUTPUT_RATE, OUTPUT_BITS, 1);
    out.write(buffer, 10*1440);
    return 0;
}