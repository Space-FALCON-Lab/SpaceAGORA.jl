variable "aws_region" {
  type        = string
  description = "AWS region for deployment"
  default     = "us-east-1"
}

variable "artifacts_bucket_name" {
  type        = string
  description = "S3 bucket for run artifacts"
}

variable "run_table_name" {
  type        = string
  description = "DynamoDB table name for run metadata"
  default     = "spaceagora-runs"
}

variable "job_queue_name" {
  type        = string
  description = "SQS queue name for simulation jobs"
  default     = "spaceagora-run-jobs"
}

variable "ecs_cluster_name" {
  type        = string
  description = "ECS cluster name"
  default     = "spaceagora-orchestrator"
}
