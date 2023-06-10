#include "encoder.hh"
#include <cstdint>

static EncoderInterface *encoder;

int main(int argc, char **argv) {

    int rate = 8000;
    int channel = 0;
    int noise_symbols = 1;
    int carrier_freq = 1500;
    int output_bits = 16;
    char out_file[] = "out.wav";

    if (argc < 3 | argc > 8) {
        std::cerr << "usage: " << argv[0] << " MESSAGE CALLSIGN [NOISE_SYMBOLS] [CARRIER_FREQUENCY] [RATE] [BITS] [CHANNEL]" << std::endl;
        return 1;
    }
    std::string mesg(argv[1]);
    const char *call_sign = argv[2];
    if (argc > 3)
        noise_symbols = std::atoi(argv[3]);
    if (argc > 4)
        carrier_freq = std::atoi(argv[4]);
    if (argc > 5)
        rate = std::atoi(argv[5]);
    if (argc > 6)
        output_bits = std::atoi(argv[6]);
    if (argc > 7)
        channel = std::atoi(argv[7]);

    char message[192] = {0};
    for (int i = 0; i < mesg.length(); i++) {
        if (i < 192)
            message[i] = mesg[i];
    }
    int channel_count = 1;
    if (channel != 0) {
        channel_count = 2;
    }
    DSP::WriteWAV<float> output_file(out_file, rate, output_bits, channel_count);
    switch (rate) {
        case 8000:
            encoder = new(std::nothrow) Encoder<8000>(&output_file);
            break;
        case 16000:
            encoder = new(std::nothrow) Encoder<16000>(&output_file);
            break;
        case 32000:
            encoder = new(std::nothrow) Encoder<32000>(&output_file);
            break;
        case 44100:
            encoder = new(std::nothrow) Encoder<44100>(&output_file);
            break;
        case 48000:
            encoder = new(std::nothrow) Encoder<48000>(&output_file);
            break;
        default:
            std::cerr << "Unsupported sample rate." ;
            std::cerr << "Supported rates: 8000/16000/32000/44100/48000."<< std::endl;
            return 1;
    }
    encoder->configure(reinterpret_cast<const uint8_t*>(&message[0]), 
        reinterpret_cast<const int8_t*>(&call_sign[0]), 
        carrier_freq, 
        noise_symbols, 
        false);
    //not sure how often to iterate. up to testing
    for (int i = 0; i < 8 + noise_symbols; ++i) {
        encoder->produce_write(channel);
	}

    return 0;
}