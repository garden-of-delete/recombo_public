#include <gtest/gtest.h>
#include "test_data_util.h"
#include "mmc_config.h"
#include "bfacf_config.h"
#include "recombo_config.h"
#include "zanalyzer_config.h"

#include <string>
#include <vector>

using namespace std;

// --- MmcConfig defaults ---

TEST(MmcConfigTest, DefaultValues)
{
    MmcConfig config;
    EXPECT_EQ(config.q, 1);
    EXPECT_DOUBLE_EQ(config.swap_ratio, 0.8);
    EXPECT_EQ(config.swap_interval, 5);
    EXPECT_EQ(config.num_samples, 0);
    EXPECT_EQ(config.steps_between_samples, 20000);
    EXPECT_EQ(config.warmup_steps, 100000);
    EXPECT_EQ(config.num_chains, 1);
    EXPECT_EQ(config.sample_mode, 's');
    EXPECT_EQ(config.seed, 0);
    EXPECT_EQ(config.sequence_type, -1);
    EXPECT_EQ(config.recombo_type, -1);
    EXPECT_FALSE(config.suppress_output);
}

// --- MmcConfig validation ---

TEST(MmcConfigTest, ValidateRequiresInputFile)
{
    MmcConfig config;
    config.output_file = "out";
    config.num_samples = 10;
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("input_file"));
}

TEST(MmcConfigTest, ValidateRequiresHaltingCriterion)
{
    MmcConfig config;
    config.input_file = "in";
    config.output_file = "out";
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("halting"));
}

TEST(MmcConfigTest, ValidateAcceptsTimeLimitAsHaltingCriterion)
{
    MmcConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.time_limit_hours = 1;
    vector<string> errors = config.validate();
    EXPECT_TRUE(errors.empty());
}

TEST(MmcConfigTest, ValidateRejectsInvalidMode)
{
    MmcConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.num_samples = 10;
    config.sample_mode = 'x';
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("sample_mode"));
}

TEST(MmcConfigTest, ValidateMultiChainRequiresZRange)
{
    MmcConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.num_samples = 10;
    config.num_chains = 5;
    config.z_min = 0.2;
    config.z_max = 0.1; // invalid: min >= max
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
}

TEST(MmcConfigTest, ValidateFilterModeRequiresRecomboParas)
{
    MmcConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.num_samples = 10;
    config.sample_mode = 'f';
    // sequence_type and recombo_type left at -1
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("recomboParas"));
}

TEST(MmcConfigTest, ValidatePassesWithValidConfig)
{
    MmcConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.num_samples = 10;
    config.z_min = 0.2;
    config.z_max = 0.2;
    vector<string> errors = config.validate();
    EXPECT_TRUE(errors.empty());
}

// --- MmcConfig parse_recombo_paras ---

TEST(MmcConfigTest, ParseRecomboParasDirect)
{
    MmcConfig config;
    EXPECT_TRUE(config.parse_recombo_paras("ds"));
    EXPECT_EQ(config.sequence_type, 1);
    EXPECT_EQ(config.recombo_type, 0);
}

TEST(MmcConfigTest, ParseRecomboParasInverted)
{
    MmcConfig config;
    EXPECT_TRUE(config.parse_recombo_paras("is"));
    EXPECT_EQ(config.sequence_type, 0);
    EXPECT_EQ(config.recombo_type, 0);
}

TEST(MmcConfigTest, ParseRecomboParasMixed)
{
    MmcConfig config;
    EXPECT_TRUE(config.parse_recombo_paras("dsp"));
    EXPECT_EQ(config.sequence_type, 1);
    EXPECT_EQ(config.recombo_type, 4);
}

TEST(MmcConfigTest, ParseRecomboParasAllCodes)
{
    struct { const char* code; int seq; int rec; } cases[] = {
        {"is", 0, 0}, {"ip", 0, 1}, {"in", 0, 2}, {"ia", 0, 3},
        {"isp", 0, 4}, {"isn", 0, 5}, {"isa", 0, 6},
        {"ds", 1, 0}, {"dp", 1, 1}, {"dn", 1, 2}, {"da", 1, 3},
        {"dsp", 1, 4}, {"dsn", 1, 5}, {"dsa", 1, 6},
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        MmcConfig config;
        EXPECT_TRUE(config.parse_recombo_paras(cases[i].code))
            << "failed to parse: " << cases[i].code;
        EXPECT_EQ(config.sequence_type, cases[i].seq) << "code: " << cases[i].code;
        EXPECT_EQ(config.recombo_type, cases[i].rec) << "code: " << cases[i].code;
    }
}

TEST(MmcConfigTest, ParseRecomboParasRejectsInvalid)
{
    MmcConfig config;
    EXPECT_FALSE(config.parse_recombo_paras("x"));
    EXPECT_FALSE(config.parse_recombo_paras("xyz"));
    EXPECT_FALSE(config.parse_recombo_paras("dx"));
    EXPECT_FALSE(config.parse_recombo_paras(""));
    EXPECT_FALSE(config.parse_recombo_paras("dsss"));
}

// --- MmcConfig from_json ---

TEST(MmcConfigTest, FromJsonAllFields)
{
    vector<string> errors;
    MmcConfig config = MmcConfig::from_json(
        string(TEST_DATA_DIR) + "/valid_mmc_config.json", errors);
    EXPECT_TRUE(errors.empty()) << errors[0];

    EXPECT_EQ(config.input_file, "initial/3_1");
    EXPECT_EQ(config.output_file, "results/trefoil");
    EXPECT_DOUBLE_EQ(config.z_min, 0.2);
    EXPECT_DOUBLE_EQ(config.z_max, 0.2);
    EXPECT_EQ(config.q, 1);
    EXPECT_DOUBLE_EQ(config.swap_ratio, 0.8);
    EXPECT_EQ(config.swap_interval, 5);
    EXPECT_EQ(config.num_samples, 5000);
    EXPECT_EQ(config.steps_between_samples, 20000);
    EXPECT_EQ(config.warmup_steps, 100000);
    EXPECT_EQ(config.num_chains, 1);
    EXPECT_EQ(config.sample_mode, 's');
    EXPECT_EQ(config.seed, 42);
    EXPECT_EQ(config.block_file_size, 1000);
    EXPECT_EQ(config.time_limit_hours, 0);
    EXPECT_TRUE(config.suppress_output);
}

TEST(MmcConfigTest, FromJsonPartialFieldsKeepDefaults)
{
    vector<string> errors;
    MmcConfig config = MmcConfig::from_json(
        string(TEST_DATA_DIR) + "/partial_mmc_config.json", errors);
    EXPECT_TRUE(errors.empty()) << errors[0];

    EXPECT_EQ(config.input_file, "initial/0_1");
    EXPECT_EQ(config.num_samples, 100);
    // Unspecified fields should keep defaults
    EXPECT_EQ(config.q, 1);
    EXPECT_DOUBLE_EQ(config.swap_ratio, 0.8);
    EXPECT_EQ(config.sample_mode, 's');
    EXPECT_EQ(config.seed, 0);
}

TEST(MmcConfigTest, FromJsonInvalidProducesErrors)
{
    vector<string> errors;
    MmcConfig config = MmcConfig::from_json(
        string(TEST_DATA_DIR) + "/invalid_mmc_config.json", errors);
    // Parse should succeed, but validate should fail
    EXPECT_TRUE(errors.empty());
    vector<string> validation = config.validate();
    EXPECT_FALSE(validation.empty());
}

TEST(MmcConfigTest, FromJsonMissingFileProducesError)
{
    vector<string> errors;
    MmcConfig::from_json("/nonexistent/path.json", errors);
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("cannot open"));
}

TEST(MmcConfigTest, FromJsonStringRecomboTypes)
{
    vector<string> errors;
    MmcConfig config = MmcConfig::from_json(
        string(TEST_DATA_DIR) + "/recombo_paras_config.json", errors);
    EXPECT_TRUE(errors.empty()) << errors[0];
    EXPECT_EQ(config.sequence_type, 1); // "direct"
    EXPECT_EQ(config.recombo_type, 0);  // "standard"
}

// --- MmcConfig to_json round-trip ---

TEST(MmcConfigTest, ToJsonRoundTrip)
{
    MmcConfig original;
    original.input_file = "initial/3_1";
    original.output_file = "results/test";
    original.z_min = 0.18;
    original.z_max = 0.21;
    original.q = 2;
    original.num_samples = 500;
    original.seed = 42;
    original.sample_mode = 'b';

    string json_str = original.to_json();

    // Parse back
    string err;
    json11::Json parsed = json11::Json::parse(json_str, err);
    EXPECT_TRUE(err.empty()) << err;
    EXPECT_EQ(parsed["input_file"].string_value(), "initial/3_1");
    EXPECT_EQ(parsed["output_file"].string_value(), "results/test");
    EXPECT_DOUBLE_EQ(parsed["zmin"].number_value(), 0.18);
    EXPECT_DOUBLE_EQ(parsed["zmax"].number_value(), 0.21);
    EXPECT_EQ(parsed["q"].int_value(), 2);
    EXPECT_EQ(parsed["n_samples"].int_value(), 500);
    EXPECT_EQ(parsed["seed"].int_value(), 42);
    EXPECT_EQ(parsed["sample_mode"].string_value(), "b");
}

// --- BfacfConfig ---

TEST(BfacfConfigTest, ValidateRequiresZAndSteps)
{
    BfacfConfig config;
    config.input_file = "in";
    config.output_file = "out";
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
}

TEST(BfacfConfigTest, ValidatePassesWithRequiredFields)
{
    BfacfConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.z = 0.2;
    config.steps = 1000;
    vector<string> errors = config.validate();
    EXPECT_TRUE(errors.empty());
}

TEST(BfacfConfigTest, ValidateRejectsInvalidFormat)
{
    BfacfConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.z = 0.2;
    config.steps = 1000;
    config.output_format = "invalid";
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("output_format"));
}

// --- RecomboConfig ---

TEST(RecomboConfigTest, ValidateRequiresArcAndParas)
{
    RecomboConfig config;
    config.input_file = "in";
    config.output_file = "out";
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
    EXPECT_GE(errors.size(), 3u); // minarc, maxarc, recomboParas
}

TEST(RecomboConfigTest, ValidatePassesWithRequiredFields)
{
    RecomboConfig config;
    config.input_file = "in";
    config.output_file = "out";
    config.min_arc = 4;
    config.max_arc = 50;
    config.parse_recombo_paras("ds");
    vector<string> errors = config.validate();
    EXPECT_TRUE(errors.empty());
}

TEST(RecomboConfigTest, ParseRecomboParasMatchesMmcConfig)
{
    // Verify both config types parse the same codes identically
    MmcConfig mmc;
    RecomboConfig rec;
    const char* codes[] = {"is", "ds", "dp", "dsp", "isa"};
    for (size_t i = 0; i < sizeof(codes)/sizeof(codes[0]); i++) {
        EXPECT_TRUE(mmc.parse_recombo_paras(codes[i]));
        EXPECT_TRUE(rec.parse_recombo_paras(codes[i]));
        EXPECT_EQ(mmc.sequence_type, rec.sequence_type) << "code: " << codes[i];
        EXPECT_EQ(mmc.recombo_type, rec.recombo_type) << "code: " << codes[i];
    }
}

// --- ZAnalyzerConfig ---

TEST(ZAnalyzerConfigTest, ValidateRejectsZminGeZmax)
{
    ZAnalyzerConfig config;
    config.input_file = "in";
    config.z_min = 0.3;
    config.z_max = 0.1;
    vector<string> errors = config.validate();
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("zmin"));
}

TEST(ZAnalyzerConfigTest, ValidatePassesWithValidRange)
{
    ZAnalyzerConfig config;
    config.input_file = "in";
    config.z_min = 0.1;
    config.z_max = 0.3;
    vector<string> errors = config.validate();
    EXPECT_TRUE(errors.empty());
}

// --- MmcConfig from_cli ---

TEST(MmcConfigTest, FromCliAllFlags)
{
    const char* argv[] = {"mmc", "initial/3_1", "results/out",
        "-zmin", "0.18", "-zmax", "0.21", "-q", "2", "-sr", "0.9",
        "-s", "10", "-n", "5000", "-c", "30000", "-m", "3",
        "-w", "500000", "-mode", "b", "-seed", "42",
        "-minarc", "4", "-maxarc", "50", "-targetlength", "120",
        "-bfs", "1000", "-t", "8", "+s", "-recomboParas", "ds"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    MmcConfig config = MmcConfig::from_cli(argc, argv, errors);
    EXPECT_TRUE(errors.empty()) << errors[0];

    EXPECT_EQ(config.input_file, "initial/3_1");
    EXPECT_EQ(config.output_file, "results/out");
    EXPECT_DOUBLE_EQ(config.z_min, 0.18);
    EXPECT_DOUBLE_EQ(config.z_max, 0.21);
    EXPECT_EQ(config.q, 2);
    EXPECT_DOUBLE_EQ(config.swap_ratio, 0.9);
    EXPECT_EQ(config.swap_interval, 10);
    EXPECT_EQ(config.num_samples, 5000);
    EXPECT_EQ(config.steps_between_samples, 30000);
    EXPECT_EQ(config.num_chains, 3);
    EXPECT_EQ(config.warmup_steps, 500000);
    EXPECT_EQ(config.sample_mode, 'b');
    EXPECT_EQ(config.seed, 42);
    EXPECT_EQ(config.min_arc, 4);
    EXPECT_EQ(config.max_arc, 50);
    EXPECT_EQ(config.target_recombo_length, 120);
    EXPECT_EQ(config.block_file_size, 1000);
    EXPECT_EQ(config.time_limit_hours, 8);
    EXPECT_TRUE(config.suppress_output);
    EXPECT_EQ(config.sequence_type, 1);
    EXPECT_EQ(config.recombo_type, 0);
}

TEST(MmcConfigTest, FromCliMinimalArgs)
{
    const char* argv[] = {"mmc", "initial/0_1", "results/out", "-n", "10"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    MmcConfig config = MmcConfig::from_cli(argc, argv, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_EQ(config.input_file, "initial/0_1");
    EXPECT_EQ(config.num_samples, 10);
    // Defaults preserved
    EXPECT_EQ(config.q, 1);
    EXPECT_EQ(config.sample_mode, 's');
}

TEST(MmcConfigTest, FromCliUnrecognizedOption)
{
    const char* argv[] = {"mmc", "in", "out", "-n", "10", "-bogus", "val"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    MmcConfig::from_cli(argc, argv, errors);
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("bogus"));
}

TEST(MmcConfigTest, FromCliRejectsConfigFlag)
{
    const char* argv[] = {"mmc", "in", "out", "-config", "file.json"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    MmcConfig::from_cli(argc, argv, errors);
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("mutually exclusive"));
}

TEST(MmcConfigTest, FromCliBadRecomboParas)
{
    const char* argv[] = {"mmc", "in", "out", "-n", "10", "-recomboParas", "xyz"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    MmcConfig::from_cli(argc, argv, errors);
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(string::npos, errors[0].find("recomboParas"));
}

// --- BfacfConfig from_cli ---

TEST(BfacfConfigTest, FromCliAllFlags)
{
    const char* argv[] = {"bfacf", "initial/0_1", "output",
        "-z", "0.2", "-q", "2", "-n", "1000", "-k", "100",
        "-bfs", "500", "-seed", "42", "-fmt", "text", "+s"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    BfacfConfig config = BfacfConfig::from_cli(argc, argv, errors);
    EXPECT_TRUE(errors.empty()) << errors[0];

    EXPECT_EQ(config.input_file, "initial/0_1");
    EXPECT_EQ(config.output_file, "output");
    EXPECT_DOUBLE_EQ(config.z, 0.2);
    EXPECT_EQ(config.q, 2);
    EXPECT_EQ(config.steps, 1000);
    EXPECT_EQ(config.save_interval, 100);
    EXPECT_EQ(config.block_file_size, 500);
    EXPECT_EQ(config.seed, 42);
    EXPECT_EQ(config.output_format, "text");
    EXPECT_TRUE(config.suppress_output);
}

TEST(BfacfConfigTest, FromCliUnrecognizedOption)
{
    const char* argv[] = {"bfacf", "in", "out", "-z", "0.2", "-n", "100", "--verbose"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    BfacfConfig::from_cli(argc, argv, errors);
    EXPECT_FALSE(errors.empty());
}

// --- RecomboConfig from_cli ---

TEST(RecomboConfigTest, FromCliAllFlags)
{
    const char* argv[] = {"recombo", "initial/3_1", "output",
        "-minarc", "4", "-maxarc", "50", "-targetlength", "120",
        "-recomboParas", "is", "-seed", "42", "-fmt", "newsud", "+s"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    RecomboConfig config = RecomboConfig::from_cli(argc, argv, errors);
    EXPECT_TRUE(errors.empty()) << errors[0];

    EXPECT_EQ(config.input_file, "initial/3_1");
    EXPECT_EQ(config.output_file, "output");
    EXPECT_EQ(config.min_arc, 4);
    EXPECT_EQ(config.max_arc, 50);
    EXPECT_EQ(config.target_length, 120);
    EXPECT_EQ(config.sequence_type, 0);
    EXPECT_EQ(config.recombo_type, 0);
    EXPECT_EQ(config.seed, 42);
    EXPECT_EQ(config.output_format, "newsud");
    EXPECT_TRUE(config.suppress_output);
}

// --- ZAnalyzerConfig from_cli ---

TEST(ZAnalyzerConfigTest, FromCliAllFlags)
{
    const char* argv[] = {"zAnalyzer", "initial/0_1",
        "-l", "60", "-tol", "5", "-zmin", "0.1", "-zmax", "0.3",
        "-q", "1", "-c", "1000", "-w", "5000", "-s", "42"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    vector<string> errors;
    ZAnalyzerConfig config = ZAnalyzerConfig::from_cli(argc, argv, errors);
    EXPECT_TRUE(errors.empty()) << errors[0];

    EXPECT_EQ(config.input_file, "initial/0_1");
    EXPECT_EQ(config.target_length, 60);
    EXPECT_EQ(config.tolerance, 5);
    EXPECT_DOUBLE_EQ(config.z_min, 0.1);
    EXPECT_DOUBLE_EQ(config.z_max, 0.3);
    EXPECT_EQ(config.q, 1);
    EXPECT_EQ(config.steps_between_samples, 1000);
    EXPECT_EQ(config.warmup, 5000);
    EXPECT_EQ(config.seed, 42);
}
