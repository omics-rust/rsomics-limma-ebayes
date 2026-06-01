use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_limma_ebayes::{Options, run, write_results};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-limma-ebayes", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// log-expression matrix TSV: header = sample ids, col 1 = gene ids.
    pub expr: PathBuf,
    /// Design matrix TSV: header = coefficient names, col 1 = sample ids.
    #[arg(long)]
    design: PathBuf,
    /// Contrast matrix TSV (contrasts.fit): col 1 = coefficient names.
    #[arg(long)]
    contrast: Option<PathBuf>,
    /// 1-based coefficient (or contrast) to tabulate; default = last.
    #[arg(long)]
    coef: Option<usize>,
    /// Expected proportion of differentially expressed genes (B-statistic).
    #[arg(long, default_value_t = 0.01)]
    proportion: f64,
    /// Results TSV destination; "-" is stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let opts = Options {
            expr: &self.expr,
            design: &self.design,
            contrast: self.contrast.as_deref(),
            coef: self.coef,
            proportion: self.proportion,
        };
        let res = run(&opts)?;

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        write_results(&res, &mut out)?;

        if !self.common.quiet {
            eprintln!(
                "{} genes, coef '{}', df.prior={:.4} df.total={:.4} s2.prior={:.6}",
                res.rows.len(),
                res.coef_name,
                res.df_prior,
                res.df_total,
                res.s2_prior
            );
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Per-gene linear model + empirical-Bayes moderated t-statistics.",
    origin: Some(Origin {
        upstream: "limma lmFit + eBayes + topTable",
        upstream_license: "GPL (>=2)",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.2202/1544-6115.1027"),
    }),
    usage_lines: &["<expr.tsv> --design <design.tsv> [--contrast <c.tsv>] [--coef N] [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "design",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Design matrix TSV (header = coefficient names, col 1 = sample ids).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "contrast",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Contrast matrix TSV; applies contrasts.fit before moderation.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "coef",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: None,
                description: "1-based coefficient/contrast to tabulate.",
                why_default: Some(
                    "Last coefficient — the typical treatment effect in an intercept+group design.",
                ),
            },
            FlagSpec {
                short: None,
                long: "proportion",
                aliases: &[],
                value: Some("<p>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.01"),
                description: "Expected proportion of DE genes (B-statistic prior).",
                why_default: Some("limma's default."),
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Results TSV destination; \"-\" is stdout.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Two-group test, tabulate the group coefficient",
            command: "rsomics-limma-ebayes E.tsv --design design.tsv --coef 2 -o top.tsv",
        },
        Example {
            description: "Contrast of two groups",
            command: "rsomics-limma-ebayes E.tsv --design design.tsv --contrast c.tsv > top.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
