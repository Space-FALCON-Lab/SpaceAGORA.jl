terraform {
  required_version = ">= 1.6.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = ">= 5.0"
    }
  }
}

provider "aws" {
  region = var.aws_region
}

resource "aws_s3_bucket" "artifacts" {
  bucket = var.artifacts_bucket_name
}

resource "aws_dynamodb_table" "runs" {
  name         = var.run_table_name
  billing_mode = "PAY_PER_REQUEST"
  hash_key     = "job_id"

  attribute {
    name = "job_id"
    type = "S"
  }
}

resource "aws_sqs_queue" "jobs" {
  name                       = var.job_queue_name
  visibility_timeout_seconds = 900
  message_retention_seconds  = 86400
}

resource "aws_cloudwatch_log_group" "orchestrator" {
  name              = "/spaceagora/orchestrator"
  retention_in_days = 14
}

resource "aws_ecs_cluster" "main" {
  name = var.ecs_cluster_name
}

resource "aws_apigatewayv2_api" "http" {
  name          = "spaceagora-orchestrator-api"
  protocol_type = "HTTP"
}

resource "aws_apigatewayv2_stage" "default" {
  api_id      = aws_apigatewayv2_api.http.id
  name        = "$default"
  auto_deploy = true
}
